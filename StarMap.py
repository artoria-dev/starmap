import colorsys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation
import matplotlib.colors as mcolors


class StarMap:
    def __init__(self, latitude, longitude, epoch_date, universal_time,
                 constellation_lines=None, constellation_labels=None):
        self.latitude = latitude
        self.longitude = longitude
        self.epoch_date = epoch_date
        self.universal_time = universal_time
        self.catalog = None
        self.animation = None
        self.is_animating = False
        self.constellation_lines = constellation_lines or {}
        self.constellation_labels = constellation_labels or {}
        self.constellation_colors = self._generate_constellation_colors()

        self.fig = plt.figure(figsize=(16, 10), facecolor='#000814')
        self.ax = self.fig.add_axes([0.05, 0.2, 0.9, 0.75], facecolor='#000000')

        self._setup_controls()

    def _generate_constellation_colors(self):
        np.random.seed(42)
        return {
            const_id: colorsys.hsv_to_rgb(np.random.random(), 0.4, 0.8)
            for const_id in self.constellation_lines
        }

    def _setup_controls(self):
        # label, min, max, curr, rect
        controls_config = [
            ('observer lat', -90, 90, self.latitude, [0.15, 0.12, 0.6, 0.02]),
            ('observer lon', -180, 180, self.longitude, [0.15, 0.09, 0.6, 0.02]),
            ('temporal offset', -365, 365, 0, [0.15, 0.06, 0.6, 0.02]),
            ('sidereal time', 0, 23.99, 12, [0.15, 0.03, 0.6, 0.02])
        ]

        self.sliders = {}
        for label, vmin, vmax, vinit, rect in controls_config:
            ax_slider = self.fig.add_axes(rect, facecolor='#000000')
            slider = Slider(ax_slider, label, vmin, vmax, valinit=vinit, color='#00b4d8')
            slider.label.set_color('white')
            slider.valtext.set_color('white')
            slider.on_changed(self._update_view)

            key = label.lower().replace(' ', '_')
            self.sliders[key] = slider

        self.base_epoch = datetime.strptime(self.epoch_date, "%Y-%m-%d")

        ax_btn = self.fig.add_axes([0.8, 0.03, 0.1, 0.05], facecolor='#000000')
        self.animate_btn = Button(ax_btn, 'animate', color='#003566', hovercolor='#00b4d8')
        self.animate_btn.label.set_color('white')
        self.animate_btn.on_clicked(self._toggle_animation)

    def _toggle_animation(self, event):
        if self.is_animating:
            self.is_animating = False
            self.animate_btn.label.set_text('animate')
            if self.animation:
                self.animation.event_source.stop()
        else:
            self.is_animating = True
            self.animate_btn.label.set_text('halt')
            # need an init render
            self._render()
            self.animation = FuncAnimation(
                self.fig,
                self._advance_time,
                interval=30,
                cache_frame_data=False
            )

    def _advance_time(self, frame):
        # need to incr sideral time for animation
        current_time = self.sliders['sidereal_time'].val
        new_time = current_time + 0.08

        if new_time >= 24.0:
            new_time -= 24.0
            day_offset = self.sliders['temporal_offset'].val + 1
            self.sliders['temporal_offset'].set_val(day_offset)

        self.sliders['sidereal_time'].set_val(new_time)

    def _update_view(self, val):
        self.latitude = self.sliders['observer_lat'].val
        self.longitude = self.sliders['observer_lon'].val

        day_offset = int(self.sliders['temporal_offset'].val)
        target_date = self.base_epoch + timedelta(days=day_offset)
        self.epoch_date = target_date.strftime("%Y-%m-%d")

        hours = int(self.sliders['sidereal_time'].val)
        minutes = int((self.sliders['sidereal_time'].val - hours) * 60)
        self.universal_time = f"{hours:02d}:{minutes:02d}"

        self._render()

    def _compute_julian_date_and_lst(self):
        dt = datetime.strptime(f"{self.epoch_date} {self.universal_time}", "%Y-%m-%d %H:%M")

        # jan/feb req
        year = dt.year - 1 if dt.month <= 2 else dt.year
        month = dt.month + 12 if dt.month <= 2 else dt.month

        # julian date
        a = int(365.25 * (year + 4716))
        b = int(30.6001 * (month + 1))
        jd = a + b + dt.day + (2 - year // 100 + year // 400) - 1524.5
        jd += (dt.hour + dt.minute / 60.0) / 24.0

        days_since_j2000 = jd - 2451545.0

        # local sideral time
        gst = (18.697374558 + 24.06570982441908 * days_since_j2000) % 24
        lst = (gst + self.longitude / 15.0) % 24

        return days_since_j2000, lst

    def _compute_sun_position(self, days_since_j2000):
        # get suns equatorial coordinates
        # mean lon
        L = (280.46 + 0.9856474 * days_since_j2000) % 360
        # mean anomaly
        g = np.radians((357.528 + 0.9856003 * days_since_j2000) % 360)
        # ecliptic longitude
        lambda_sun = np.radians(L + 1.915 * np.sin(g) + 0.02 * np.sin(2 * g))
        # obliquity of ecliptic
        epsilon = np.radians(23.439 - 0.0000004 * days_since_j2000)

        # convert to equatorial coords
        ra = np.degrees(np.arctan2(np.cos(epsilon) * np.sin(lambda_sun), np.cos(lambda_sun)))
        dec = np.degrees(np.arcsin(np.sin(epsilon) * np.sin(lambda_sun)))

        return (ra / 15.0) % 24, dec

    def _compute_altitude(self, ra_hours, dec_degrees, lst):
        lat_rad = np.radians(self.latitude)
        dec_rad = np.radians(dec_degrees)
        ha_rad = np.radians((lst - ra_hours) * 15)

        sin_alt = np.sin(lat_rad) * np.sin(dec_rad) + \
                  np.cos(lat_rad) * np.cos(dec_rad) * np.cos(ha_rad)

        return np.arcsin(sin_alt)

    def _render(self):
        if self.catalog is None:
            return

        days, lst = self._compute_julian_date_and_lst()
        sun_ra, sun_dec = self._compute_sun_position(days)

        self.ax.clear()
        self._draw_twilight(sun_ra, sun_dec)
        star_colors = self._compute_star_colors(lst)
        star_sizes = self._compute_star_sizes(lst)
        self._draw_stars(star_colors, star_sizes)
        self._draw_constellations(lst)
        self._draw_horizon(lst)
        self.ax.axvline(lst, color='#ff4d6d', ls='--', lw=0.8, alpha=0.3, label='local meridian')
        # sun
        self.ax.scatter(sun_ra, sun_dec, s=200, c='#ffc300', edgecolors='#ffb703',
                        marker='o', zorder=5, label='solar position')
        self._configure_plot()

    def _draw_twilight(self, sun_ra, sun_dec):
        # fancy atmospheric twilight gradient
        ra_grid, dec_grid = np.meshgrid(
            np.linspace(0, 24, 72),
            np.linspace(-90, 90, 36)
        )
        # angular distance from sun
        cos_dist = np.sin(np.radians(dec_grid)) * np.sin(np.radians(sun_dec)) + \
                   np.cos(np.radians(dec_grid)) * np.cos(np.radians(sun_dec)) * \
                   np.cos(np.radians((ra_grid - sun_ra) * 15))

        angular_dist = np.degrees(np.arccos(np.clip(cos_dist, -1, 1)))
        twilight = np.clip((108 - angular_dist) / 18, 0, 1)

        cmap = mcolors.LinearSegmentedColormap.from_list(
            'dawn', [(0, 0, 0, 0), (0.1, 0.2, 0.5, 0.2)]
        )
        self.ax.contourf(ra_grid, dec_grid, twilight, levels=12, cmap=cmap, zorder=1)

    def _compute_star_colors(self, lst):
        altitudes = np.degrees(self._compute_altitude(
            self.catalog['ra_hours'].values,
            self.catalog['dec'].values,
            lst
        ))

        colors = []
        white = np.array([1.0, 1.0, 1.0])
        sunset = np.array([1.0, 0.4, 0.1])
        underground = np.array([0.2, 0.3, 0.4])

        for alt in altitudes:
            if alt >= 15:
                colors.append((*white, 0.9))
            elif alt >= 0:
                blend = alt / 15.0
                color = sunset + (white - sunset) * blend
                colors.append((*color, 0.9))
            else:
                alpha = np.clip(0.4 + (alt / 90.0) * 0.3, 0.1, 0.4)
                colors.append((*underground, alpha))

        return colors

    def _compute_star_sizes(self, lst):
        altitudes = self._compute_altitude(
            self.catalog['ra_hours'].values,
            self.catalog['dec'].values,
            lst
        )

        alt_deg = np.degrees(altitudes)
        visible = alt_deg > 0

        # stmospheric extinction
        airmass = np.ones_like(alt_deg)
        airmass[visible] = 1.0 / np.sin(np.maximum(altitudes[visible], np.radians(1.5)))

        apparent_mag = self.catalog['phot_g_mean_mag'].values.copy()
        apparent_mag[visible] += airmass[visible] * 0.25
        apparent_mag[~visible] += 2.5

        return np.clip(14 * (8.5 - apparent_mag) / (apparent_mag + 1), 0.1, 14)

    def _draw_stars(self, colors, sizes):
        # outer glow
        self.ax.scatter(
            self.catalog['ra_hours'],
            self.catalog['dec'],
            s=sizes * 5,
            c=colors,
            alpha=0.04,
            edgecolors='none',
            zorder=2
        )
        # core
        self.ax.scatter(
            self.catalog['ra_hours'],
            self.catalog['dec'],
            s=sizes,
            c=colors,
            edgecolors='none',
            zorder=3,
            label='stellar objects'
        )

    def _draw_constellations(self, lst):
        for const_id, paths in self.constellation_lines.items():
            color = self.constellation_colors.get(const_id, (0.5, 0.5, 0.5))

            for path in paths:
                ra_pts = np.array([p[0] for p in path])
                dec_pts = np.array([p[1] for p in path])
                altitudes = np.degrees(self._compute_altitude(ra_pts, dec_pts, lst))

                # split at RA discontinuities ?
                splits = np.where(np.abs(np.diff(ra_pts)) > 12)[0] + 1

                for segment in np.split(np.arange(len(ra_pts)), splits):
                    if len(segment) == 0:
                        continue

                    visible = altitudes[segment] > 0
                    alpha = 0.2 if np.any(visible) else 0.08
                    linewidth = 0.6 if np.any(visible) else 0.4

                    self.ax.plot(ra_pts[segment], dec_pts[segment],
                                 color=color, lw=linewidth, alpha=alpha)

        # labels
        for const_id, meta in self.constellation_labels.items():
            altitude = np.degrees(self._compute_altitude(meta['ra'], meta['dec'], lst))
            visible = altitude > 0

            alpha = 0.5 if visible else 0.15
            color = self.constellation_colors.get(const_id, 'white') if visible else '#457b9d'

            # always show major, minor only when visible
            if meta['rank'] <= (1 if not visible else 2):
                text = meta['name'] if meta['rank'] == 1 else meta['abbr']
                weight = 'bold' if meta['rank'] == 1 else 'normal'

                self.ax.text(meta['ra'], meta['dec'], text,
                             color=color, fontsize=7, alpha=alpha,
                             ha='center', va='center', fontweight=weight)

    def _draw_horizon(self, lst):
        dec_range = np.linspace(-90, 90, 400)

        # hour angle at horizon
        cos_ha = -np.tan(np.radians(self.latitude)) * np.tan(np.radians(dec_range))
        valid = np.abs(cos_ha) <= 1

        ha = np.degrees(np.arccos(cos_ha[valid])) / 15.0

        # conv to RA
        ra_west = (lst - ha) % 24
        ra_east = (lst + ha) % 24

        self.ax.scatter(
            np.concatenate([ra_west, ra_east]),
            np.concatenate([dec_range[valid], dec_range[valid]]),
            s=0.8, c='#ffd60a', alpha=0.3, label='local horizon'
        )

    def _configure_plot(self):
        self.ax.set_title(
            f"star map | epoch: {self.epoch_date} | utc: {self.universal_time}",
            color='white', fontsize=12, pad=10
        )
        self.ax.set(xlim=(24, 0), ylim=(-90, 90),
                    xlabel="RA", ylabel="declination")

        self.ax.set_xticks(np.arange(0, 25, 2))
        self.ax.set_xticklabels([f"{int(x)}h" for x in np.arange(0, 25, 2)])

        self.ax.xaxis.label.set_color('white')
        self.ax.yaxis.label.set_color('white')
        self.ax.tick_params(colors='white', labelsize=8)

        self.ax.legend(loc='upper right', facecolor='#000814',
                       edgecolor='white', labelcolor='white',
                       fontsize=7, framealpha=0.4)

        self.fig.canvas.draw_idle()

    def show(self, catalog_df=None):
        self.catalog = catalog_df if catalog_df is not None else pd.read_feather('catalog.feather')
        self._render()
        plt.show()