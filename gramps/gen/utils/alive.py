#
# Gramps - a GTK+/GNOME based genealogy program
#
# Copyright (C) 2000-2007  Donald N. Allingham
# Copyright (C) 2009       Gary Burton
# Copyright (C) 2011       Tim G L Lyons
# Copyright (C) 2024       Cameron Davidson
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""
A utility to make a best guess if a person is alive.  This is used to provide
privacy in reports and exports.
"""

# -------------------------------------------------------------------------
#
# Standard python modules
#
# -------------------------------------------------------------------------
import logging

# -------------------------------------------------------------------------
#
# Gramps modules
#
# -------------------------------------------------------------------------
from ..display.name import displayer as name_displayer
from ..lib.date import Date, Today
from ..lib.person import Person
from ..errors import DatabaseError
from ..const import GRAMPS_LOCALE as glocale
from ..proxy.proxybase import ProxyDbBase

LOG = logging.getLogger(".gen.utils.alive")

_ = glocale.translation.sgettext

# -------------------------------------------------------------------------
#
# Constants from config .ini keys
#
# -------------------------------------------------------------------------
# cache values; use refresh_constants() if they change
try:
    from ..config import config

    _MAX_AGE_PROB_ALIVE = config.get("behavior.max-age-prob-alive")
    _MAX_SIB_AGE_DIFF = config.get("behavior.max-sib-age-diff")
    _AVG_GENERATION_GAP = config.get("behavior.avg-generation-gap")
    _MIN_GENERATION_YEARS = config.get("behavior.min-generation-years")
except ImportError:
    # Utils used as module not part of GRAMPS
    _MAX_AGE_PROB_ALIVE = 110
    _MAX_SIB_AGE_DIFF = 20
    _AVG_GENERATION_GAP = 20
    _MIN_GENERATION_YEARS = 13


# -------------------------------------------------------------------------
#
# ProbablyAlive class
#
# -------------------------------------------------------------------------
class ProbablyAlive:
    """
    An object to hold the parameters for considering someone alive.
    """

    def __init__(
        self,
        db,
        max_sib_age_diff=None,
        max_age_prob_alive=None,
        avg_generation_gap=None,
        min_generation_years=None,
    ):
        self.db = db
        if max_sib_age_diff is None:
            max_sib_age_diff = _MAX_SIB_AGE_DIFF
        if max_age_prob_alive is None:
            max_age_prob_alive = _MAX_AGE_PROB_ALIVE
        if avg_generation_gap is None:
            avg_generation_gap = _AVG_GENERATION_GAP
        if min_generation_years is None:
            min_generation_years = _MIN_GENERATION_YEARS
        self.MAX_SIB_AGE_DIFF = max_sib_age_diff
        self.MAX_AGE_PROB_ALIVE = max_age_prob_alive
        self.AVG_GENERATION_GAP = avg_generation_gap
        self.MIN_GENERATION_YEARS = min_generation_years
        self.pset = set()

    def probably_alive_range(self, person, is_spouse=False):
        """
        Find likely birth and death date ranges, either from dates of actual
        events recorded in the db or else estimating range outer limits from
        other events in their lives or those of close family.

        Returns: (birth_date, death_date, explain_text, related_person)
        """
        # where appropriate, some derived dates are expressed as a range.
        if person is None:
            return (None, None, "", None)
        self.pset = set()
        birth_date = None
        death_date = None
        known_to_be_dead = False
        min_birth_year = None
        max_birth_year = None
        min_birth_year_from_death = None  # values derived from 110 year extrapolations
        max_birth_year_from_death = None
        # these min/max parameters are simply years
        sib_birth_min, sib_birth_max = (None, None)
        explain_birth_min = ""
        explain_birth_max = ""
        explain_death = ""

        def get_person_bdm(class_or_handle):
            """
                Looks up birth and death events for referenced person,
                using fallback dates if necessary.
                The dates will always be either None or valid values, avoiding EMPTYs

            returns  (birth_date, death_date, death_found, explain_birth, explain_death)
                                 for the referenced person
            """
            birth_date = None
            death_date = None
            death_found = False
            explain_birth = ""
            explain_death = ""

            if not class_or_handle:
                return (
                    birth_date,
                    death_date,
                    death_found,
                    explain_birth,
                    explain_death,
                )

            if isinstance(class_or_handle, Person):
                thisperson = class_or_handle
            elif isinstance(class_or_handle, str):
                thisperson = self.db.get_person_from_handle(class_or_handle)
            else:
                thisperson = None

            if not thisperson:
                LOG.debug("    get_person_bdm: null person called")
                return (
                    birth_date,
                    death_date,
                    death_found,
                    explain_birth,
                    explain_death,
                )
            # is there an actual death record?  Even if yes, there may be no date,
            # in which case the EMPTY date is reported for the event.
            death_ref = thisperson.get_death_ref()
            if death_ref and death_ref.get_role().is_primary():
                evnt = self.db.get_event_from_handle(death_ref.ref)
                if evnt:
                    death_found = True
                    dateobj = evnt.get_date_object()
                    if dateobj and dateobj.is_valid():
                        death_date = dateobj
                        explain_death = _("date")

            # at this stage death_date is None or a valid date.
            # death_found is true if person is known to be dead,
            #        whether or not a date was found.
            # If we have no death_date then look for fallback even such as Burial.
            # These fallbacks are fairly good indications that someone's not alive.
            # If that date itself is not valid, it means we know they are dead
            #  but not when they died. So keep checking in case we get a date.
            if not death_date:
                for ev_ref in person.get_primary_event_ref_list():
                    if ev_ref:
                        evnt = self.db.get_event_from_handle(ev_ref.ref)
                        if evnt and evnt.type.is_death_fallback():
                            death_date_fb = evnt.get_date_object()
                            death_found = True
                            if death_date_fb.is_valid():
                                death_date = death_date_fb
                                explain_death = _("date fallback")
                                if death_date.get_modifier() == Date.MOD_NONE:
                                    death_date.set_modifier(Date.MOD_BEFORE)
                                break  # we found a valid date, stop looking.
            # At this point:
            # * death_found is False: (no death indication found); or
            # * death_found is True. (death confirmed somehow);  In which case:
            #       * (death_date is valid) some form of death date found; or
            #       * (death_date is None and no date was recorded)
            # now repeat, looking for birth date
            birth_ref = thisperson.get_birth_ref()
            if birth_ref and birth_ref.get_role().is_primary():
                evnt = self.db.get_event_from_handle(birth_ref.ref)
                if evnt:
                    dateobj = evnt.get_date_object()
                    if dateobj and dateobj.is_valid():
                        birth_date = dateobj
                        explain_birth = _("date")

            # to here:
            #   birth_date is None: either no birth record or else no date reported; or
            #   birth_date is a valid date
            # Look for Baptism, etc events.
            # These are fairly good indications of someone's birth date.
            if not birth_date:
                for ev_ref in person.get_primary_event_ref_list():
                    evnt = self.db.get_event_from_handle(ev_ref.ref)
                    if evnt and evnt.type.is_birth_fallback():
                        birth_date_fb = evnt.get_date_object()
                        if birth_date_fb and birth_date_fb.is_valid():
                            birth_date = birth_date_fb
                            explain_birth = _("date fallback")
                            break

            return (birth_date, death_date, death_found, explain_birth, explain_death)

        birth_date, death_date, known_to_be_dead, explain_birth_min, explain_death = (
            get_person_bdm(person)
        )

        explanation = (
            _("DIRECT birth: ") + explain_birth_min + _(", death: ") + explain_death
        )
        if death_date is not None and birth_date is not None:
            return (birth_date, death_date, explanation, person)  # direct self evidence

        # birth and/or death dates are not known, so let's see what we can estimate.
        # First: minimum is X years before death;
        # Second: get the parent's birth/death dates if available, so we can constrain
        # to sensible values - mother's age and parent's death.
        # Finally: get birth dates for any full siblings to further constrain.
        # Currently only look at full siblings, ranges would get wider for half sibs.

        if birth_date is None:
            # only need to estimate birth_date if we have no more direct evidence.
            if death_date is not None:
                # person died so guess initial limits to birth date
                if death_date.get_year_valid():
                    max_birth_year_from_death = death_date.get_year()
                    min_birth_year_from_death = max_birth_year_from_death - self.MAX_AGE_PROB_ALIVE

            m_birth, m_death = (None, None)  # mother's birth and death dates
            f_birth, f_death = (None, None)  # father's
            parents = None  # Family with parents
            parenth = person.get_main_parents_family_handle()
            if parenth:
                parents = self.db.get_family_from_handle(parenth)
                mum = parents.get_mother_handle()
                m_birth, m_death = get_person_bdm(mum)[0:2]
                dad = parents.get_father_handle()
                f_birth, f_death = get_person_bdm(dad)[0:2]
            # now scan siblings
            family_list = person.get_parent_family_handle_list()
            for family_handle in family_list:
                family = self.db.get_family_from_handle(family_handle)
                if family is None:
                    continue
                if parents is not None and family != parents:
                    LOG.debug(
                        "      skipping family %s.",
                        family.get_gramps_id(),
                    )
                    continue
                for child_ref in family.get_child_ref_list():
                    child_handle = child_ref.ref
                    child = self.db.get_person_from_handle(child_handle)
                    if child is None or child == person:
                        continue
                    need_birth_fallback = True
                    # Go through once looking for direct evidence:
                    # extract the range of birth dates, either direct or fallback
                    for ev_ref in child.get_primary_event_ref_list():
                        evnt = self.db.get_event_from_handle(ev_ref.ref)
                        if evnt and evnt.type.is_birth():
                            dobj = evnt.get_date_object()
                            if dobj and dobj.get_year_valid():
                                year = dobj.get_year()
                                need_birth_fallback = False
                                if sib_birth_min is None or year < sib_birth_min:
                                    sib_birth_min = year
                                if sib_birth_max is None or year > sib_birth_max:
                                    sib_birth_max = year
                    # scan even list again looking for fallback:
                    if need_birth_fallback:
                        for ev_ref in child.get_primary_event_ref_list():
                            evnt = self.db.get_event_from_handle(ev_ref.ref)
                            if evnt and evnt.type.is_birth_fallback():
                                dobj = evnt.get_date_object()
                                if dobj and dobj.get_year_valid():
                                    # if sibling birth date too far away, then
                                    # cannot be alive:
                                    year = dobj.get_year()
                                    if sib_birth_min is None or year < sib_birth_min:
                                        sib_birth_min = year
                                    if sib_birth_max is None or year > sib_birth_max:
                                        sib_birth_max = year
            # Now evaluate estimate based on parents and siblings:
            # Make sure child is born after both parents are old enough
            if m_birth:
                min_birth_year = m_birth.get_year() + self.MIN_GENERATION_YEARS
                explain_birth_min = _("mother's age")
            if f_birth:
                min_from_f = f_birth.get_year() + self.MIN_GENERATION_YEARS
                if min_birth_year is None or min_from_f > min_birth_year:
                    min_birth_year = min_from_f
                    explain_birth_min = _("father's age")
            if min_birth_year_from_death:
                if min_birth_year is None or min_birth_year_from_death > min_birth_year:
                    min_birth_year = min_birth_year_from_death
                    explain_birth_min = _("from death date")
            # Calculate the latest year that the child could have been born
            if m_death:
                max_birth_year = m_death.get_year()
                explain_birth_max = _("mother's death")
            if f_death:
                max_from_f = f_death.get_year() + 1
                if max_birth_year is None or max_from_f < max_birth_year:
                    max_birth_year = max_from_f
                    explain_birth_max = _("father's death")
            if max_birth_year_from_death:
                if max_birth_year is None or max_birth_year_from_death < max_birth_year:
                    max_birth_year = max_birth_year_from_death
                    explain_birth_max = _("person's death")

        # sib_xx_min/max are either both None or both have a value (maybe the same)
        if sib_birth_max:
            min_from_sib = sib_birth_max - self.MAX_SIB_AGE_DIFF
            if min_birth_year is None or min_from_sib > min_birth_year:
                min_birth_year = min_from_sib
                explain_birth_min = _("oldest sibling's age")

            max_from_sib = sib_birth_min + self.MAX_SIB_AGE_DIFF
            if max_birth_year is None or max_from_sib < max_birth_year:
                max_birth_year = max_from_sib
                explain_birth_max = _("youngest sibling's age")

        if birth_date is None or not birth_date.is_valid():
            birth_date = Date()  # make sure we have an empty date
            #  use proxy estimate
            if min_birth_year and max_birth_year:
                # create a range set
                birth_range = list(Date.EMPTY + Date.EMPTY)
                birth_range[Date._POS_YR] = min_birth_year
                birth_range[Date._POS_RYR] = max_birth_year
                birth_date.set(modifier=Date.MOD_RANGE, value=tuple(birth_range))
            else:
                if min_birth_year:
                    birth_date.set_yr_mon_day(min_birth_year, 1, 1)
                    birth_date.set_modifier(Date.MOD_AFTER)
                elif max_birth_year:
                    birth_date.set_yr_mon_day(max_birth_year, 12, 31)
                    birth_date.set_modifier(Date.MOD_BEFORE)
            birth_date.recalc_sort_value()

        # If we have no death date but we know death has happened then
        # we set death range somewhere between birth and yesterday.
        # otherwise we assume MAX years after birth
        if death_date is None:
            if birth_date and birth_date.is_valid():
                death_date = Date(birth_date)
                max_death_date = birth_date.copy_offset_ymd(
                    year=self.MAX_AGE_PROB_ALIVE
                )
                if known_to_be_dead:
                    if max_death_date.match(Today(), ">="):
                        max_death_date = Today()
                        max_death_date.set_yr_mon_day_offset(
                            day=-1
                        )  # make it yesterday
                    # range start value stays at birth date
                    death_date.set_modifier(Date.MOD_RANGE)
                    death_date.set_text_value("")
                    death_date.set2_yr_mon_day(
                        max_death_date.get_year(),
                        max_death_date.get_month(),
                        max_death_date.get_day(),
                    )
                    explain_death = _("birth date and known to be dead")
                else:
                    death_date.set_yr_mon_day_offset(year=self.MAX_AGE_PROB_ALIVE)
                    if death_date.is_compound():
                        death_date.set2_yr_mon_day_offset(year=self.MAX_AGE_PROB_ALIVE)
                    explain_death = _("birth date")
                death_date.recalc_sort_value()
            else:
                death_date = Date()

        # at this stage we should have valid dates for both birth and death,
        #  or else both are zero (if None then it's a bug).
        if explain_birth_max == "":
            explanation = _("birth: ") + explain_birth_min
        else:
            explanation = (
                _("birth: ") + explain_birth_min + _(" and ") + explain_birth_max
            )
        explanation += _(", death: ") + explain_death
        explanation = "2ND + " + explanation
        if birth_date.is_valid() and death_date.is_valid():
            return (birth_date, death_date, explanation, person)

        # Try looking for descendants that were born more than a lifespan
        # ago.

        def descendants_too_old(person, years):
            """
            Recursively scan descendants' tree to determine likely birth/death
            dates for the person specified.
            years: gets incremented by average generation gap as we descend the tree.
            """
            # FIXME:  this code does not follow an optimum path.
            # currently: for each child of the ref. person:
            #   if it has an actual birth date - return estimated parent dates
            #   if it has an actual death date - return estimated parent dates
            #     recurse to that child's children
            #   if recursion returns empty, then
            #   check fallback death dates for child.
            # Problems - this is close to a depth-first search, while it should
            # instead check all children, then all grandchildren, etc.
            # The cost of fixing probably exceeds the value gained.

            # First, a couple of routines to reduce duplication
            def get_bd_from_child_birth_event(evnt, explain):
                """
                Estimate birth/death dates of a person from a child/grandchild's birth
                (or birth fallback) date
                Returns: (birth_date, death_date, explain, child)
                        if valid dates were found,
                        otherwise: None
                """
                nonlocal child, years
                dobj = evnt.get_date_object()
                retval = None
                if dobj.get_start_date() != Date.EMPTY:
                    dretbirth = Date(dobj)
                    dretbirth.set_yr_mon_day_offset(year=(-years))
                    dretdeath = dretbirth.copy_offset_ymd(self.MAX_AGE_PROB_ALIVE)
                    if dretbirth.is_compound():
                        # if is not very meaningful to adjust the upper limit of a
                        # compound date, however it is far worse to leave any unchanged.
                        # A better alternative might be to remove the 2nd date.
                        dretbirth.set2_yr_mon_day_offset(year=(-years))
                        dretdeath.set2_yr_mon_day_offset(
                            year=(-(self.MAX_AGE_PROB_ALIVE + years))
                        )
                    retval = (dretbirth, dretdeath, explain, child)
                return retval

            def get_bd_from_child_death_event(evnt, explain):
                """
                Estimate birth/death dates from a child's death (or fallback)
                date, having already decided there is no useful birth date.
                This process is very uncertain, as the child may have died in
                infancy or at the age of 100.
                Returns: (birth_date, death_date, explain, child)
                        if valid dates were found,
                        otherwise: None
                """
                nonlocal child, years
                dobj = evnt.get_date_object()
                if dobj.get_start_date() != Date.EMPTY:
                    dretbirth = Date(dobj)
                    child_death_yr = dobj.get_year()
                    min_birth_yr = child_death_yr - (years + self.MAX_AGE_PROB_ALIVE)
                    max_birth_yr = child_death_yr - years
                    birth_range = list(Date.EMPTY + Date.EMPTY)
                    birth_range[Date._POS_YR] = min_birth_yr
                    birth_range[Date._POS_RYR] = max_birth_yr
                    dretbirth.set(modifier=Date.MOD_RANGE, value=tuple(birth_range))
                    dretdeath = dretbirth.copy_offset_ymd(self.MAX_AGE_PROB_ALIVE)
                    return (
                        dretbirth,
                        dretdeath,
                        explain,
                        child,
                    )
                return None

            if person.handle in self.pset:
                return (None, None, "", None)
            self.pset.add(person.handle)
            for family_handle in person.get_family_handle_list():
                family = self.db.get_family_from_handle(family_handle)
                if not family:
                    # can happen with LivingProxyDb(PrivateProxyDb(db))
                    continue
                for child_ref in family.get_child_ref_list():
                    child_handle = child_ref.ref
                    child = self.db.get_person_from_handle(child_handle)
                    child_birth_ref = child.get_birth_ref()
                    if child_birth_ref:
                        child_birth = self.db.get_event_from_handle(child_birth_ref.ref)
                        rval = get_bd_from_child_birth_event(
                            child_birth,
                            _("descendant birth date"),
                        )
                        if rval is not None:
                            return rval
                    # Check birth fallback data:
                    for ev_ref in child.get_primary_event_ref_list():
                        evnt = self.db.get_event_from_handle(ev_ref.ref)
                        if evnt and evnt.type.is_birth_fallback():
                            rval = get_bd_from_child_birth_event(
                                evnt,
                                _("descendant birth-related date"),
                            )
                            if rval is not None:
                                return rval
                    # check primary death events
                    child_death_ref = child.get_death_ref()
                    if child_death_ref:
                        child_death = self.db.get_event_from_handle(child_death_ref.ref)
                        rval = get_bd_from_child_death_event(
                            child_death, _("descendant death date")
                        )
                        if rval is not None:
                            return rval
                    # Check death fallback data:
                    for ev_ref in child.get_primary_event_ref_list():
                        evnt = self.db.get_event_from_handle(ev_ref.ref)
                        if evnt and evnt.type.is_death_fallback():
                            rval = get_bd_from_child_death_event(
                                evnt, _("descendant death-related date")
                            )
                            if rval is not None:
                                return rval
                    date1, date2, explain, other = descendants_too_old(
                        child, years + self.AVG_GENERATION_GAP
                    )
                    if date1 and date2:
                        return date1, date2, explain, other

            return (None, None, "", None)

        # If there are descendants that are too old for the person to have
        # been alive in the current year then they must be dead.

        date1, date2, explain, other = None, None, "", None
        try:
            date1, date2, explain, other = descendants_too_old(
                person, self.AVG_GENERATION_GAP
            )
        except RuntimeError:
            raise DatabaseError(
                _("Database error: loop in %s's descendants")
                % name_displayer.display(person)
            )

        if date1 and date2:
            return (date1, date2, explain, other)

        def ancestors_too_old(person, year):
            if person.handle in self.pset:
                return (None, None, "", None)
            self.pset.add(person.handle)
            LOG.debug(
                "ancestors_too_old('%s', %d)", name_displayer.display(person), year
            )
            family_handle = person.get_main_parents_family_handle()
            if family_handle:
                family = self.db.get_family_from_handle(family_handle)
                if not family:
                    # can happen with LivingProxyDb(PrivateProxyDb(db))
                    return (None, None, "", None)
                father_handle = family.get_father_handle()
                if father_handle:
                    father = self.db.get_person_from_handle(father_handle)
                    father_birth_ref = father.get_birth_ref()
                    if father_birth_ref and father_birth_ref.get_role().is_primary():
                        father_birth = self.db.get_event_from_handle(
                            father_birth_ref.ref
                        )
                        dobj = father_birth.get_date_object()
                        if dobj.get_start_date() != Date.EMPTY:
                            return (
                                dobj.copy_offset_ymd(-year),
                                dobj.copy_offset_ymd(-year + self.MAX_AGE_PROB_ALIVE),
                                _("ancestor birth date"),
                                father,
                            )
                    father_death_ref = father.get_death_ref()
                    if father_death_ref and father_death_ref.get_role().is_primary():
                        father_death = self.db.get_event_from_handle(
                            father_death_ref.ref
                        )
                        dobj = father_death.get_date_object()
                        if dobj.get_start_date() != Date.EMPTY:
                            return (
                                dobj.copy_offset_ymd(-year - self.MAX_AGE_PROB_ALIVE),
                                dobj.copy_offset_ymd(
                                    -year
                                    - self.MAX_AGE_PROB_ALIVE
                                    + self.MAX_AGE_PROB_ALIVE
                                ),
                                _("ancestor death date"),
                                father,
                            )

                    # Check fallback data:
                    for ev_ref in father.get_primary_event_ref_list():
                        evnt = self.db.get_event_from_handle(ev_ref.ref)
                        if evnt and evnt.type.is_birth_fallback():
                            dobj = evnt.get_date_object()
                            if dobj.get_start_date() != Date.EMPTY:
                                return (
                                    dobj.copy_offset_ymd(-year),
                                    dobj.copy_offset_ymd(
                                        -year + self.MAX_AGE_PROB_ALIVE
                                    ),
                                    _("ancestor birth-related date"),
                                    father,
                                )

                        elif evnt and evnt.type.is_death_fallback():
                            dobj = evnt.get_date_object()
                            if dobj.get_start_date() != Date.EMPTY:
                                return (
                                    dobj.copy_offset_ymd(
                                        -year - self.MAX_AGE_PROB_ALIVE
                                    ),
                                    dobj.copy_offset_ymd(
                                        -year
                                        - self.MAX_AGE_PROB_ALIVE
                                        + self.MAX_AGE_PROB_ALIVE
                                    ),
                                    _("ancestor death-related date"),
                                    father,
                                )

                    date1, date2, explain, other = ancestors_too_old(
                        father, year - self.AVG_GENERATION_GAP
                    )
                    if date1 and date2:
                        return date1, date2, explain, other

                mother_handle = family.get_mother_handle()
                if mother_handle:
                    mother = self.db.get_person_from_handle(mother_handle)
                    mother_birth_ref = mother.get_birth_ref()
                    if mother_birth_ref and mother_birth_ref.get_role().is_primary():
                        mother_birth = self.db.get_event_from_handle(
                            mother_birth_ref.ref
                        )
                        dobj = mother_birth.get_date_object()
                        if dobj.get_start_date() != Date.EMPTY:
                            return (
                                dobj.copy_offset_ymd(-year),
                                dobj.copy_offset_ymd(-year + self.MAX_AGE_PROB_ALIVE),
                                _("ancestor birth date"),
                                mother,
                            )
                    mother_death_ref = mother.get_death_ref()
                    if mother_death_ref and mother_death_ref.get_role().is_primary():
                        mother_death = self.db.get_event_from_handle(
                            mother_death_ref.ref
                        )
                        dobj = mother_death.get_date_object()
                        if dobj.get_start_date() != Date.EMPTY:
                            return (
                                dobj.copy_offset_ymd(-year - self.MAX_AGE_PROB_ALIVE),
                                dobj.copy_offset_ymd(
                                    -year
                                    - self.MAX_AGE_PROB_ALIVE
                                    + self.MAX_AGE_PROB_ALIVE
                                ),
                                _("ancestor death date"),
                                mother,
                            )

                    # Check fallback data:
                    for ev_ref in mother.get_primary_event_ref_list():
                        evnt = self.db.get_event_from_handle(ev_ref.ref)
                        if evnt and evnt.type.is_birth_fallback():
                            dobj = evnt.get_date_object()
                            if dobj.get_start_date() != Date.EMPTY:
                                return (
                                    dobj.copy_offset_ymd(-year),
                                    dobj.copy_offset_ymd(
                                        -year + self.MAX_AGE_PROB_ALIVE
                                    ),
                                    _("ancestor birth-related date"),
                                    mother,
                                )

                        elif evnt and evnt.type.is_death_fallback():
                            dobj = evnt.get_date_object()
                            if dobj.get_start_date() != Date.EMPTY:
                                return (
                                    dobj.copy_offset_ymd(
                                        -year - self.MAX_AGE_PROB_ALIVE
                                    ),
                                    dobj.copy_offset_ymd(
                                        -year
                                        - self.MAX_AGE_PROB_ALIVE
                                        + self.MAX_AGE_PROB_ALIVE
                                    ),
                                    _("ancestor death-related date"),
                                    mother,
                                )

                    date1, date2, explain, other = ancestors_too_old(
                        mother, year - self.AVG_GENERATION_GAP
                    )
                    if date1 and date2:
                        return (date1, date2, explain, other)

            return (None, None, "", None)

        try:
            # If there are ancestors that would be too old in the current year
            # then assume our person must be dead too.
            date1, date2, explain, other = ancestors_too_old(
                person, -int(self.AVG_GENERATION_GAP)
            )
        except RuntimeError:
            raise DatabaseError(
                _("Database error: loop in %s's ancestors")
                % name_displayer.display(person)
            )
        if date1 and date2:
            return (date1, date2, explain, other)

        # final test is against spouse details, which involves a recursive call
        # to probably_alive_range().
        # This test gives higher uncertainty than others ...
        # We allow for an age difference +/- AVG_GENERATION_GAP
        # which, assuming defaults, results in 150 year "probably alive" range.
        # In reality,

        if not is_spouse:  # if you are not in recursion, let's recurse:
            LOG.debug(
                "    ----- trying spouse check: birth %s, death %s",
                "valid" if birth_date.is_valid() else "indeterminate",
                "valid" if death_date.is_valid() else "indeterminate",
            )
            for family_handle in person.get_family_handle_list():
                family = self.db.get_family_from_handle(family_handle)
                if family:
                    mother_handle = family.get_mother_handle()
                    father_handle = family.get_father_handle()
                    spouse = None
                    if mother_handle == person.handle and father_handle:
                        spouse = self.db.get_person_from_handle(father_handle)
                    elif father_handle == person.handle and mother_handle:
                        spouse = self.db.get_person_from_handle(mother_handle)
                    if spouse is not None:
                        date1, date2, explain, other = self.probably_alive_range(
                            spouse, is_spouse=True
                        )
                        if date1 and date1.get_year() != 0:
                            return (
                                Date().copy_ymd(
                                    date1.get_year() - self.AVG_GENERATION_GAP
                                ),
                                Date().copy_ymd(
                                    date1.get_year()
                                    + self.AVG_GENERATION_GAP
                                    + self.MAX_AGE_PROB_ALIVE
                                ),
                                _("a spouse's birth-related date, ") + explain,
                                other,
                            )
                        if date2 and date2.get_year() != 0:
                            return (
                                Date().copy_ymd(
                                    date2.get_year()
                                    - self.AVG_GENERATION_GAP
                                    - self.MAX_AGE_PROB_ALIVE
                                ),
                                Date().copy_ymd(
                                    date2.get_year() + self.AVG_GENERATION_GAP
                                ),
                                _("a spouse's death-related date, ") + explain,
                                other,
                            )

                    # Let's check the family events and see if we find something
                    for ref in family.get_event_ref_list():
                        if ref:
                            event = self.db.get_event_from_handle(ref.ref)
                            if event:
                                date = event.get_date_object()
                                year = date.get_year()
                                if year != 0:
                                    other = None
                                    if person.handle == mother_handle and father_handle:
                                        other = self.db.get_person_from_handle(
                                            father_handle
                                        )
                                    elif (
                                        person.handle == father_handle and mother_handle
                                    ):
                                        other = self.db.get_person_from_handle(
                                            mother_handle
                                        )
                                    return (
                                        Date().copy_ymd(year - self.AVG_GENERATION_GAP),
                                        Date().copy_ymd(
                                            year
                                            - self.AVG_GENERATION_GAP
                                            + self.MAX_AGE_PROB_ALIVE
                                        ),
                                        _("event with spouse"),
                                        other,
                                    )

        # If we can't find any reason to believe that they are dead we
        # must assume they are alive.

        return (None, None, "", None)


# -------------------------------------------------------------------------
#
# probably_alive
#
# -------------------------------------------------------------------------
def probably_alive(
    person,
    db,
    current_date=None,
    limit=0,
    max_sib_age_diff=None,
    max_age_prob_alive=None,
    avg_generation_gap=None,
    return_range=False,
):
    """
    Return true if the person may be alive on current_date.

    This works by a process of elimination. If we can't find a good
    reason to believe that someone is dead then we assume they must
    be alive.

    :param current_date: a date object that is not estimated or modified
                         (defaults to today)
    :param limit: number of years to check beyond death_date
    :param max_sib_age_diff: maximum sibling age difference, in years
    :param max_age_prob_alive: maximum age of a person, in years
    :param avg_generation_gap: average generation gap, in years
    """
    LOG.debug(
        " === [%s] %s: ",
        person.get_gramps_id(),
        person.get_primary_name().get_gedcom_name(),
    )
    # First, get the probable birth and death ranges for
    # this person from the real database:
    birth, death, explain, relative = probably_alive_range(
        person, db, max_sib_age_diff, max_age_prob_alive, avg_generation_gap
    )
    if current_date is None:
        current_date = Today()
    elif not current_date.is_valid():
        current_date = Today()

    if not explain.startswith("DIRECT"):
        if relative is  None:
            rel_id = "nobody"
        else:
            rel_id = relative.get_gramps_id()
        LOG.debug(
            "      b.%s, d.%s vs %s - %s to [%s]",
            birth,
            death,
            current_date,
            explain,
            rel_id,
        )
    if not birth or not death:
        # no evidence, must consider alive
        LOG.debug(
            "      [%s] %s: decided alive - no evidence",
            person.get_gramps_id(),
            person.get_primary_name().get_gedcom_name(),
        )
        return (True, None, None, _("no evidence"), None) if return_range else True
    # must have dates from here:
    if limit:
        death += limit  # add these years to death
    # Finally, check to see if current_date is between dates
    # ---true if  current_date >= birth(min)   and  true if current_date < death
    # these include true if current_date is within the estimated range
    result = current_date.match(birth, ">=") and current_date.match(death, "<")
    if not explain.startswith("DIRECT"):
        (bthmin, bthmax) = birth.get_start_stop_range()
        (dthmin, dthmax) = death.get_start_stop_range()
        (dmin, dmax) = current_date.get_start_stop_range()
        LOG.debug(
            "        alive=%s, btest: %s, dtest: %s (born %s-%s, dd %s-%s) vs (%s-%s)",
            result,
            current_date.match(birth, ">="),
            current_date.match(death, "<"),
            bthmin,
            bthmax,
            dthmin,
            dthmax,
            dmin,
            dmax,
        )
    if return_range:
        return (result, birth, death, explain, relative)

    return result


def probably_alive_range(
    person, db, max_sib_age_diff=None, max_age_prob_alive=None, avg_generation_gap=None
):
    """
    Computes estimated birth and death date ranges.
    Returns: (birth_date, death_date, explain_text, related_person)
    """
    # First, find the real database to use all people
    # for determining alive status:
    basedb = db
    while isinstance(basedb, ProxyDbBase):
        basedb = basedb.db
    # Now, we create a wrapper for doing work:
    pbac = ProbablyAlive(
        basedb, max_sib_age_diff, max_age_prob_alive, avg_generation_gap
    )
    return pbac.probably_alive_range(person)


def update_constants():
    """
    Used to update the constants that are cached in this module.
    """

    global _MAX_AGE_PROB_ALIVE, _MAX_SIB_AGE_DIFF
    global _AVG_GENERATION_GAP, _MIN_GENERATION_YEARS
    _MAX_AGE_PROB_ALIVE = config.get("behavior.max-age-prob-alive")
    _MAX_SIB_AGE_DIFF = config.get("behavior.max-sib-age-diff")
    _AVG_GENERATION_GAP = config.get("behavior.avg-generation-gap")
    _MIN_GENERATION_YEARS = config.get("behavior.min-generation-years")
