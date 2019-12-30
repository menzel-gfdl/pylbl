from collections import namedtuple
from html.parser import HTMLParser

from numpy import float32


Isotopologue = namedtuple("Isotopologue", ["id", "abundance", "mass"])


class HitranIsotopologueHTMLParser(HTMLParser):
    """Parses the isotopologue properties for a molecule from the HITRAN website.

    Attributes:
        abundance: The abundance of the current isotopologue.
        column: Index of the current column in the HTML table.
        formula: Chemical formula of the current molecule.
        id: HITRAN global id for the isotopologue.
        inside_formula: Flag telling if a molecule is currently being processed.
        inside_row: Flag telling if the parser is currently inside a table row.
        iso_list: List of Isotopologue namedtuples for the current molecule.
        isotopologues: Dictionary mapping a molecule to its list of isotopologue namedtuples.
        mass: Mass of the current isotopologue.
    """
    def __init__(self):
        HTMLParser.__init__(self)
        self.abundance = ""
        self.column = 0
        self.formula = ""
        self.id = ""
        self.inside_formula = False
        self.inside_row = False
        self.iso_list = []
        self.isotopologues = {}
        self.mass = ""

    def handle_data(self, data):
        """Stores data based on the parser's current state.

        Args:
            data: The data contained in the current HTML tag.
        """
        if self.inside_formula:
            self.formula += data.strip()
        elif self.inside_row:
            if self.column == 1:
                self.id += data.strip()
            elif self.column == 5:
                self.abundance += data.strip()
            elif self.column == 6:
                self.mass += data.strip()

    def handle_endtag(self, tag):
        """Resets flags and stores data based on the parser's current state.

        Args:
            tag: Current HTML tag.
        """
        if tag == "h4" and self.inside_formula:
            self.inside_formula = False
            self.formula = self.formula.split(":")[-1].strip()
        elif tag == "tr" and self.inside_row and self.column > 0:
            abundance = self.abundance.replace("\xa0{}\xa010".format(chr(215)), "e")
            self.abundance = ""
            self.iso_list.append(Isotopologue(int(self.id), abundance, float32(self.mass)))
            self.id = ""
            self.mass = ""
            self.inside_row = False
        elif tag == "tbody":
            self.isotopologues[self.formula] = self.iso_list
            self.formula = ""
            self.iso_list = []

    def handle_starttag(self, tag, attrs):
        """Sets flags based on the parser's current state.

        Args:
            tag: Current HTML tag.
            attrs: List of tuples containing HTML tag attributes.
        """
        if tag == "h4":
            self.inside_formula = True
        elif tag == "tr":
            self.column = 0
            self.inside_row = True
        elif tag == "td" and self.inside_row:
            self.column += 1


class HitranMoleculeIdHTMLParser(HTMLParser):
    """Parses HITRAN molecule ids from the HITRAN website.

    Attributes:
        formula: Chemical formula of the current molecule.
        id: HITRAN integer identifier of the current molecule.
        inside_row: Flag telling if the parser is inside a row of a HTML table.
        inside_table: Flag telling if the parser is inside the correct HTML table.
        molecule_ids: Dictionary mapping chemical formulae to their HITRAN ids.
    """
    def __init__(self):
        HTMLParser.__init__(self)
        self.formula = ""
        self.id = ""
        self.inside_row = False
        self.inside_table = False
        self.molecule_ids = {}

    def handle_data(self, data):
        """Stores the current HITRAN id or chemical symbol, based on the parser's current state.

        Args:
            data: The data contained in the current HTML tag.
        """
        if self.inside_row:
            if self.column == 1:
                self.id += data.strip()
            elif self.column == 2:
                self.formula += data.strip()

    def handle_endtag(self, tag):
        """Resets flags and stores data based on the parser's current state.

        Args:
           tag: Current HTML tag.
        """
        if tag == "table" and self.inside_table:
            self.inside_table = False
        if tag == "tr" and self.inside_row and self.column > 0:
            self.molecule_ids[self.formula] = int(self.id)
            self.formula = ""
            self.id = ""
            self.inside_row = False

    def handle_starttag(self, tag, attrs):
        """Sets flags based on the parser's current state.

        Args:
            tag: Current HTML tag.
            attrs: List of tuples containing HTML tag attributes.
        """
        if tag == "table" and ("class", "list-table") in attrs:
            self.inside_table = True
        elif tag == "tr" and self.inside_table:
            self.column = 0
            self.inside_row = True
        elif tag == "td" and self.inside_row:
            self.column += 1
