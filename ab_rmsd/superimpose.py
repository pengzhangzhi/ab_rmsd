import torch
from torch import Tensor
from typing import List


class KabschRMSD:
    @torch.no_grad()
    def rmsd(self, coords_pred: Tensor, coords: Tensor) -> Tensor:
        """
        Calculate the RMSD between two coordinates (without superimpose).
        Args:
            coords_pred: (...,3)
            coords: (...,3)
        Returns:
            rmsd: (1)
        """
        rmsd = torch.sqrt(torch.mean(torch.sum(((coords_pred - coords) ** 2), dim=-1)))
        return rmsd

    def __call__(self, coords_pred: Tensor, coords: Tensor, superimpose=False) -> Tensor:
        """
        Calculate the KabschRMSD between two coordinates.

        Args:
            coords_pred: (N, 3)
            coords: (N, 3)
        Returns:
            rmsd: (1)
        """

        if superimpose:
            rot, tran = self.calc_superimpose_transformation(coords, coords_pred)
            coords_pred = self.apply_transformation(rot, tran, coords_pred)
        
        return self.rmsd(coords, coords_pred)

    @torch.no_grad()
    def apply_transformation(self, rotation, translation, coords):
        """
        Apply transformation to coordinates.
        Args:
            rotation: (3,3)
            translation: (1,3)
            coords: (N,3)
        Returns:
            transformed coordinates: (N,3)
        """
        assert rotation.shape == (
            3,3,
        ), "rotation matrix should be 3x3 but got {}".format(rotation.shape)
        assert translation.shape == (
            1,3,
        ), "translation matrix should be 1x3 but got {}".format(translation.shape)
        assert (
            len(coords.shape) == 2 and coords.shape[1] == 3
        ), "coords should be Nx3 but got {}".format(coords.shape)
        return (rotation @ coords.t()).t() + translation

    @torch.no_grad()
    def calc_superimpose_transformation(self, coords_tgt, coords_src):
        """
        Calculate the transformation from source coordinates to target coordinates.

        Args:
            coords_tgt: target coordinates, [N, 3]
            coords_src: source coordinates, [N, 3]
        Returns:
            rotation: (3,3)
            translation: (1,3)
        """
        assert (
            len(coords_tgt.shape) == 2 and coords_tgt.shape[1] == 3
        ), "coords_tgt should be Nx3 but got {}".format(coords_tgt.shape)
        assert (
            len(coords_src.shape) == 2 and coords_src.shape[1] == 3
        ), "coords_src should be Nx3 but got {}".format(coords_src.shape)

        coords_pred_mean = coords_tgt.mean(dim=0, keepdim=True)  # (1,3)
        coords_mean = coords_src.mean(dim=0, keepdim=True)  # (1,3)

        A = (coords_tgt - coords_pred_mean).transpose(0, 1) @ (coords_src - coords_mean)

        U, S, Vt = torch.linalg.svd(A)

        corr_mat = torch.diag(
            torch.tensor([1, 1, torch.sign(torch.det(A))], device=coords_tgt.device)
        )
        rotation = (U @ corr_mat) @ Vt
        translation = coords_pred_mean - torch.t(rotation @ coords_mean.t())  # (1,3)

        return rotation, translation
