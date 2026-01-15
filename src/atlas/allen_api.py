"""Allen Brain Atlas SDK wrapper for CCF integration."""

import numpy as np
from pathlib import Path

# Allen SDK imports
try:
    from allensdk.core.reference_space_cache import ReferenceSpaceCache
    ALLEN_SDK_AVAILABLE = True
except ImportError:
    ALLEN_SDK_AVAILABLE = False
    print("Warning: allensdk not installed. Some features will be limited.")


# Cache directory for Allen SDK data
CACHE_DIR = Path(__file__).parent.parent.parent / "data" / "allen_cache"


class AllenAtlas:
    """Wrapper for Allen Brain Atlas data access."""

    def __init__(self, resolution=25):
        """
        Initialize Allen Atlas wrapper.

        Args:
            resolution: CCF resolution in microns (10, 25, 50, or 100)
        """
        self.resolution = resolution
        self._annotation = None
        self._template = None
        self._structure_tree = None

        if ALLEN_SDK_AVAILABLE:
            CACHE_DIR.mkdir(parents=True, exist_ok=True)
            self.reference_space_cache = ReferenceSpaceCache(
                resolution=resolution,
                cache=True,
                manifest=str(CACHE_DIR / 'manifest.json')
            )
        else:
            self.reference_space_cache = None

    @property
    def annotation(self):
        """Get the CCF annotation volume (structure IDs)."""
        if self._annotation is None and self.reference_space_cache:
            self._annotation, _ = self.reference_space_cache.get_annotation_volume()
        return self._annotation

    @property
    def template(self):
        """Get the CCF template volume (average brain)."""
        if self._template is None and self.reference_space_cache:
            self._template, _ = self.reference_space_cache.get_template_volume()
        return self._template

    @property
    def structure_tree(self):
        """Get the structure tree for ID lookups."""
        if self._structure_tree is None and self.reference_space_cache:
            self._structure_tree = self.reference_space_cache.get_structure_tree()
        return self._structure_tree

    def get_structure_mask(self, structure_id):
        """
        Get binary mask for a brain structure.

        Args:
            structure_id: Allen CCF structure ID

        Returns:
            3D numpy array (binary mask)
        """
        if self.annotation is None:
            return None

        # Get all descendant structure IDs (structure + children)
        if self.structure_tree:
            descendant_ids = self.structure_tree.descendant_ids([structure_id])[0]
        else:
            descendant_ids = [structure_id]

        # Create mask including all descendants
        mask = np.isin(self.annotation, descendant_ids)
        return mask

    def get_structure_coords(self, structure_id):
        """
        Get all voxel coordinates for a structure.

        Args:
            structure_id: Allen CCF structure ID

        Returns:
            Nx3 array of (z, y, x) coordinates in voxels
        """
        mask = self.get_structure_mask(structure_id)
        if mask is None:
            return None

        coords = np.array(np.where(mask)).T  # Shape: (N, 3)
        return coords

    def get_structure_centroid(self, structure_id):
        """Get the centroid of a structure in voxel coordinates."""
        coords = self.get_structure_coords(structure_id)
        if coords is None or len(coords) == 0:
            return None
        return coords.mean(axis=0)

    def voxel_to_um(self, coords):
        """Convert voxel coordinates to micrometers."""
        return coords * self.resolution

    def get_structure_info(self, structure_id):
        """Get structure metadata."""
        if self.structure_tree is None:
            return None
        try:
            return self.structure_tree.get_structures_by_id([structure_id])[0]
        except (IndexError, KeyError):
            return None

    def search_structures(self, query):
        """Search for structures by name."""
        if self.structure_tree is None:
            return []
        all_structures = self.structure_tree.get_structures_by_set_id([167587189])
        matches = []
        for s in all_structures:
            if query.lower() in s['name'].lower() or query.lower() in s['acronym'].lower():
                matches.append({
                    'id': s['id'],
                    'name': s['name'],
                    'acronym': s['acronym']
                })
        return matches


# Singleton instance
_atlas = None


def get_atlas(resolution=25):
    """Get or create the singleton Allen Atlas instance."""
    global _atlas
    if _atlas is None:
        _atlas = AllenAtlas(resolution=resolution)
    return _atlas
