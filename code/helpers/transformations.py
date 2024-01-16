from typing import Iterator
import numpy as np
from ..classes.aminoacid import Aminoacid
from ..classes.protein import Protein


def check_rotation(protein_1: Protein, protein_2: Protein, depth: int) -> bool:
    """
    Check if two sets of coordinates from Protein instances are rotations of each other.

    Parameters:
    - protein_1 (Protein): The first Protein instance.
    - protein_2 (Protein): The second Protein instance.
    - depth (int): The depth or the number of Aminoacid objects to consider in the check.

    Returns:
    - bool: True if the sets are rotations of each other, False otherwise.
    """
    centered_set1, _ = _get_set_and_centroid(protein_1, depth)
    centered_set2, centroid2 = _get_set_and_centroid(protein_2, depth)

    covariance_matrix = centered_set1.T @ centered_set2

    U, _, Vt = np.linalg.svd(covariance_matrix)
    rotation_matrix = U @ Vt

    if not _is_valid_rotation_matrix(rotation_matrix):
        return False

    rotated_set1 = centered_set1 @ rotation_matrix + centroid2

    return np.allclose(rotated_set1, centered_set2)


def check_mirror(protein_1: Protein, protein_2: Protein, depth: int) -> bool:
    """
    Check if two sets of coordinates from Protein instances are mirrors of each other.

    Parameters:
    - protein_1 (Protein): The first Protein instance.
    - protein_2 (Protein): The second Protein instance.
    - depth (int): The depth or the number of Aminoacid objects to consider in the check.

    Returns:
    - bool: True if the sets are mirrors of each other, False otherwise.
    """
    centered_set1, _ = _get_set_and_centroid(protein_1, depth)
    centered_set2, centroid2 = _get_set_and_centroid(protein_2, depth)

    mirror_matrix = centered_set1.T @ centered_set2

    if not _is_valid_mirror_matrix(mirror_matrix):
        return False

    mirrored_set1 = centered_set1 @ mirror_matrix + centroid2

    return np.allclose(mirrored_set1, centered_set2)


def _traverse_linked_list(start_aminoacid: Aminoacid, depth: int) -> Iterator[Aminoacid]:
    """
    Traverse a linked list up to a specified depth.

    Parameters:
    - start_aminoacid (Aminoacid): The starting Aminoacid in the linked list.
    - depth (int): The depth or the number of Aminoacid objects to traverse.

    Yields:
    - Aminoacid: A generator yielding Aminoacid objects from the linked list.
    """
    current = start_aminoacid
    for _ in range(depth):
        if current is None:
            break
        yield current
        current = current.link


def _get_set_and_centroid(protein: Protein, depth: int):
    """
    Helper function to get coordinates set and centroid of a Protein up to a specified depth.

    Parameters:
    - protein (Protein): The Protein instance.
    - depth (int): The depth or the number of Aminoacid objects to consider.

    Returns:
    - Tuple[np.ndarray, np.ndarray]: The centered set of coordinates and the centroid.
    """
    aminoacids = _traverse_linked_list(protein.get_list(), depth)
    coordinates = np.array([aminoacid.position for aminoacid in aminoacids])
    centroid = np.mean(coordinates, axis=0)
    centered_set = coordinates - centroid
    return centered_set, centroid


def _is_valid_rotation_matrix(matrix: np.ndarray) -> bool:
    """
    Check if a matrix is a valid rotation matrix.

    Parameters:
    - matrix (np.ndarray): The matrix to check.

    Returns:
    - bool: True if the matrix is a valid rotation matrix, False otherwise.
    """
    return np.allclose(matrix @ matrix.T, np.eye(3)) and np.allclose(np.linalg.det(matrix), 1.0)


def _is_valid_mirror_matrix(matrix: np.ndarray) -> bool:
    """
    Check if a matrix is a valid mirror matrix.

    Parameters:
    - matrix (np.ndarray): The matrix to check.

    Returns:
    - bool: True if the matrix is a valid mirror matrix, False otherwise.
    """
    return np.allclose(matrix @ matrix.T, np.eye(3)) and np.allclose(np.linalg.det(matrix), -1.0)
