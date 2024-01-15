from typing import Optional, Tuple


class Aminoacid:
    """
    Represents an amino acid in a protein structure.

    Attributes
    ----------
    position : Tuple[int, int, int]
        The position of the amino acid in 3D space.
    predecessor : Optional[Aminoacid]
        The predecessor amino acid in the protein sequence. Defaults to None.
    link : Optional[Aminoacid]
        The successor amino acid in the protein sequence. Defaults to None.
    _type : str
        The type of the amino acid.

    Methods
    -------
    __init__(type: str, predecessor: Optional[Aminoacid] = None, link: Optional[Aminoacid] = None)
        Initialize an Aminoacid object with the given type and optional predecessor and link.

    get_type() -> str
        Get the type of the amino acid.

    stability_score(other: Aminoacid) -> int
        Calculate the stability score of the amino acid pair.

    __str__() -> str
        Return the string representation of the amino acid.
    """

    def __init__(self, type: str, predecessor: Optional["Aminoacid"] = None,
                 link: Optional["Aminoacid"] = None) -> None:
        """
        Initialize an Aminoacid object with the given type and optional predecessor and link.

        Parameters
        ----------
        type : str
            The type of the amino acid.
        predecessor : Optional[Aminoacid], optional
            The predecessor amino acid in the protein sequence. Defaults to None.
        link : Optional[Aminoacid], optional
            The successor amino acid in the protein sequence. Defaults to None.
        """
        self.position: Tuple[int, int, int] = (0, 0, 0)
        self.predecessor: Optional["Aminoacid"] = predecessor
        self.link: Optional["Aminoacid"] = link
        self._type: str = type

    def get_type(self) -> str:
        """
        Get the type of the amino acid.

        Returns
        -------
        str
            The type of the amino acid.
        """
        return self._type

    def stability_score(self, other: "Aminoacid") -> int:
        """
        Calculate the stability score of the amino acid pair.

        Parameters
        ----------
        other : Aminoacid
            The other amino acid to calculate the stability score with.

        Returns
        -------
        int
            The calculated stability score.
        """
        if self.get_type() == "P" or other.get_type() == "P":
            return 0

        if self.get_type() == "H" or other.get_type() == "H":
            return -1

        if self.get_type() == "C":
            return -5

    def __str__(self) -> str:
        """
        Return the string representation of the amino acid.

        Returns
        -------
        str
            The type of the amino acid.
        """
        return self._type
