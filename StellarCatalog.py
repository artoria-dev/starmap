from pathlib import Path
import pandas as pd
from astroquery.gaia import Gaia


class StellarCatalog:
    def __init__(self, limit=5000, magnitude_threshold=7.5, cache_path='catalog.feather'):
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        Gaia.ROW_LIMIT = -1
        self.limit = limit
        self.magnitude_threshold = magnitude_threshold
        self.cache_path = Path(cache_path)

    def download(self):
        print(f"fetching {self.limit} stars (mag < {self.magnitude_threshold})...")

        query = f"""
            SELECT TOP {self.limit} ra, dec, phot_g_mean_mag
            FROM gaiadr3.gaia_source
            WHERE phot_g_mean_mag < {self.magnitude_threshold}
            ORDER BY phot_g_mean_mag ASC
        """

        job = Gaia.launch_job(query)
        df = job.get_results().to_pandas()
        df['ra_hours'] = df['ra'] / 15.0

        df.to_feather(self.cache_path)
        print(f"saved {len(df)} stars to {self.cache_path}")

        return df

    def load(self):
        if self.cache_path.exists():
            print(f"loading cat from {self.cache_path}")
            return pd.read_feather(self.cache_path)

        print("no cache, downloading...")
        return self.download()