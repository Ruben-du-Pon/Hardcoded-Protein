from typing import Iterator
import numpy as np
from ..classes.aminoacid import Aminoacid
from ..classes.protein import Protein


def traverse_linked_list(start_aminoacid: Aminoacid, depth: int) -> Iterator[Aminoacid]:
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
    def is_valid_rotation_matrix(matrix: np.ndarray) -> bool:
        """
        Check if a matrix is a valid rotation matrix.

        Parameters:
        - matrix (np.ndarray): The matrix to check.

        Returns:
        - bool: True if the matrix is a valid rotation matrix, False otherwise.
        """
        return np.allclose(np.dot(matrix, matrix.T), np.eye(3)) and np.allclose(np.linalg.det(matrix), 1.0)

    # Extract coordinates from the linked list of Aminoacid objects
    set1 = np.array([aminoacid.position for aminoacid in traverse_linked_list(
        next(protein_1.get_list()), depth)])
    set2 = np.array([aminoacid.position for aminoacid in traverse_linked_list(
        next(protein_2.get_list()), depth)])

    # Find the centroid of each set
    centroid1 = np.mean(set1, axis=0)
    centroid2 = np.mean(set2, axis=0)

    # Center the sets at the origin
    centered_set1 = set1 - centroid1
    centered_set2 = set2 - centroid2

    # Compute the covariance matrix
    covariance_matrix = np.dot(centered_set1.T, centered_set2)

    # Use Singular Value Decomposition (SVD) to find the rotation matrix
    U, _, Vt = np.linalg.svd(covariance_matrix)
    rotation_matrix = np.dot(U, Vt)

    # Check if the matrix is a valid rotation matrix
    if not is_valid_rotation_matrix(rotation_matrix):
        return False

    # Apply the rotation matrix to set1
    rotated_set1 = np.dot(centered_set1, rotation_matrix) + centroid2

    # Check if rotated_set1 matches set2
    return np.allclose(rotated_set1, set2)


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
    def is_valid_mirror_matrix(matrix: np.ndarray) -> bool:
        """
        Check if a matrix is a valid mirror matrix.

        Parameters:
        - matrix (np.ndarray): The matrix to check.

        Returns:
        - bool: True if the matrix is a valid mirror matrix, False otherwise.
        """
        return np.allclose(np.dot(matrix, matrix.T), np.eye(3)) and np.allclose(np.linalg.det(matrix), -1.0)

    # Extract coordinates from the linked list of Aminoacid objects
    set1 = np.array([aminoacid.position for aminoacid in traverse_linked_list(
        next(protein_1.get_list()), depth)])
    set2 = np.array([aminoacid.position for aminoacid in traverse_linked_list(
        next(protein_2.get_list()), depth)])

    # Find the centroid of each set
    centroid1 = np.mean(set1, axis=0)
    centroid2 = np.mean(set2, axis=0)

    # Center the sets at the origin
    centered_set1 = set1 - centroid1
    centered_set2 = set2 - centroid2

    # Compute the mirror matrix
    mirror_matrix = np.dot(centered_set1.T, centered_set2)

    # Check if the matrix is a valid mirror matrix
    if not is_valid_mirror_matrix(mirror_matrix):
        return False

    # Apply the mirror matrix to set1
    mirrored_set1 = np.dot(centered_set1, mirror_matrix) + centroid2

    # Check if mirrored_set1 matches set2
    return np.allclose(mirrored_set1, set2)
