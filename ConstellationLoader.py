import json
import urllib.request
import numpy as np
import os


class ConstellationLoader:
    LINES_URL = "https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/constellations.lines.json"
    NAMES_URL = "https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/constellations.json"

    def __init__(self, cache_dir='constellation_data'):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)

    def _download_json(self, url, filename):
        cache_path = os.path.join(self.cache_dir, filename)

        if os.path.exists(cache_path):
            print(f"{filename} exists, loading from cache")
            with open(cache_path, 'r') as f:
                return json.load(f)

        print(f"downloading {filename}...")
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())

        with open(cache_path, 'w') as f:
            json.dump(data, f)

        return data

    def _convert_coords(self, lon, lat):
        if lon >= 0:
            ra = lon / 15.0
        else:
            ra = (360 + lon) / 15.0

        dec = lat

        return ra, dec

    def load_lines(self):
        data = self._download_json(self.LINES_URL, 'constellations.lines.json')

        lines_by_constellation = {}

        for feature in data['features']:
            const_id = feature['id']
            multilinestring = feature['geometry']['coordinates']
            const_lines = []
            for linestring in multilinestring:
                ra_dec_pairs = []
                for lon, lat in linestring:
                    ra, dec = self._convert_coords(lon, lat)
                    ra_dec_pairs.append((ra, dec))
                const_lines.append(ra_dec_pairs)

            lines_by_constellation[const_id] = const_lines

        return lines_by_constellation

    def load_names(self):
        data = self._download_json(self.NAMES_URL, 'constellations.json')

        labels = {}

        for feature in data['features']:
            const_id = feature['id']
            properties = feature['properties']
            coords = feature['geometry']['coordinates']

            ra, dec = self._convert_coords(coords[0], coords[1])

            name = properties.get('en', properties.get('name', const_id))

            rank = int(properties.get('rank', 2))

            labels[const_id] = {
                'name': name,
                'abbr': const_id,
                'ra': ra,
                'dec': dec,
                'rank': rank
            }

        return labels

    def load(self):
        lines = self.load_lines()
        labels = self.load_names()
        return lines, labels

    def get_major_constellations(self, labels, max_rank=1):
        return {k: v for k, v in labels.items() if v['rank'] <= max_rank}