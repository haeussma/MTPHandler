import itertools as it

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from MTPHandler.dataclasses import Plate


def visualize_plate(
    plate: Plate,
    zoom: bool = False,
    wavelengths: list[float] = [],
    static: bool = False,
):
    """Visualize a plate with all its wells and measurements."""

    if zoom:
        shared_yaxes = False
    else:
        shared_yaxes = True

    if not wavelengths:
        wavelengths = [plate.wells[0].measurements[0].wavelength]

    if not isinstance(wavelengths, list):
        wavelengths = [wavelengths]

    fig = make_subplots(
        rows=8,
        cols=12,
        shared_xaxes=True,
        subplot_titles=_generate_possible_well_ids(),
        shared_yaxes=shared_yaxes,
    )
    colors = px.colors.qualitative.Plotly

    for well in plate.wells:
        for measurement, color in zip(well.measurements, colors):
            if measurement.wavelength not in wavelengths:
                continue

            fig.add_trace(
                go.Scatter(
                    x=plate.times,
                    y=measurement.absorption,
                    name=f"{measurement.wavelength} nm",
                    mode="lines",
                    showlegend=False,
                    line=dict(color=color),
                    hovertemplate="%{y:.2f}<br>",
                ),
                col=well.x_pos + 1,
                row=well.y_pos + 1,
            )

    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)

    fig.update_layout(
        plot_bgcolor="white",
        hovermode="x",
        title=dict(
            text=f"{plate.temperatures[0]} {plate.temperature_unit.name}",
        ),
        margin=dict(l=20, r=20, t=100, b=20),
    )

    if static:
        fig.show("png")

    fig.show()


def _generate_possible_well_ids() -> list[str]:
    characters = "ABCDEFGH"
    integers = range(1, 13)  # 1 to 12

    sub_char = characters[:8]
    sub_int = integers[:12]

    # Generate combinations of characters and integers
    combinations = ["".join(item) for item in it.product(sub_char, map(str, sub_int))]

    return combinations
