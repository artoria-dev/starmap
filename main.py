from datetime import datetime
from zoneinfo import ZoneInfo

from StellarCatalog import StellarCatalog
from ConstellationLoader import ConstellationLoader
from StarMap import StarMap


def main():
    # load catalog
    catalog = StellarCatalog(limit=6000, magnitude_threshold=8.0)
    star_data = catalog.load()

    # get constellations
    loader = ConstellationLoader()
    constellation_lines, constellation_labels = loader.load()


    berlin_time = datetime.now(ZoneInfo('Europe/Berlin'))
    projection = StarMap(
        latitude=52.52,
        longitude=13.405,
        epoch_date=berlin_time.strftime('%Y-%m-%d'),
        universal_time=berlin_time.strftime('%H:%M'),
        constellation_lines=constellation_lines,
        constellation_labels=constellation_labels
    )
    projection.show(star_data)


if __name__ == '__main__':
    main()