# The sketchpad: draw on the map and see it compute

The sketchpad is a short menu in the Explorer's map tools (the pencil under the brush). Each
tool draws something on the globe, runs a real computation in your browser on it, and draws
the answer back on the map: a wave travelling down the rivers, a flattened flood, a field of
contours. A small card holds the numbers, and on a phone it is a sheet at the bottom. Every
drawing and every run is an entry in the action log, with its undo. Nothing needs a key or a
server: the package runs in the page's worker (Pyodide), and the same functions are a CLI
verb and an MCP tool.

Part of [#664](https://github.com/Rekin226/aquascope/issues/664), "Make the map compute".
The modules are experimental (see [stability.md](stability.md)).

## Draw a storm

Draw an area, say how much rain falls on it and over how many hours, and watch the runoff
travel down every river below it ([#670](https://github.com/Rekin226/aquascope/issues/670)).

![The sketchpad's menu](img/sketchpad/storm-menu-desktop.jpg)

**Draw a storm** in the menu starts the drawing: tap the area's corners and press **Done**
(or Enter, or a double tap), or drag a box. The drawing tool is the one levees are drawn with
on Flood it's valley ([floodsim.md](floodsim.md#levees-and-breaches)), with an area instead of
a line.

![The area drawn over mid-Wales, with the rain, the hours and the runoff coefficient](img/sketchpad/storm-area-desktop.jpg)

The card then asks for three numbers:

- **Rain**, in mm over the area (above 0, at most 2,000).
- **Over**, the hours it falls over (1 to 72 here; the package takes up to 720).
- **Runs off**, the runoff coefficient: the share of the rain that runs off, 0.3 by default.
  It is said plainly on the card because it is plain: a single number, the same on every
  catchment, with no soil, no wetness before the storm and no snow behind it.

**Route the storm** runs `aquascope.routing.route_storm` in the worker, the same function as
`aquascope route storm` and the MCP tool `route_storm`:

1. Every GEOGLOWS v2 reach whose unit catchment's centre (GEOGLOWS's ERA5 weight table) falls
   inside the area gets the runoff coefficient times the rain over its own unit catchment,
   spread evenly over the hours.
2. That runoff is routed down every river below with Muskingum, GEOGLOWS's own k and x per
   reach, hour by hour, until the wave has passed the last reach (Route it's engine; see
   [Route a flood](explorer.md#route-a-flood) and its check against GEOGLOWS's own
   retrospective).

![The storm's wave 16 hours after the rain started, on its way down the rivers below the area](img/sketchpad/storm-wave-desktop.jpg)

The wave plays on the rivers with the time bar's own clock (hours since the rain started),
coloured and widened by the share of the biggest peak the storm adds that is there at that
hour, through Route it's own layers. **Back to dates** shows every reach at its own peak. A
click on a reach opens its card: the wave's front, its peak and the hour of the peak.

The card says, all from the package:

- the biggest flow the storm adds anywhere and the hour it peaks;
- the hydrographs of the reach that gets the most and of the main outlet, hour by hour;
- how many reaches and km² are under the storm, the runoff in m³, and the main outlet's peak;
- **Left by the outlets**: the share of the runoff that has left by the outlets by the end of
  the hours routed, the run's own mass balance;
- the rule used, in words ("0.3 of the rain runs off: 15 mm over each unit catchment under the
  storm, spread evenly over 6 h. A plain coefficient, no soil, wetness or snow.").

**Drag a corner** of the area and the storm is routed again from where it is now. The move
and the new run are each in the action log; undo takes back the run, then the move.

![A corner dragged east: 114 reaches under the storm now, routed again](img/sketchpad/storm-moved-desktop.jpg)

![The action log: the area, the run, the move and the run again, each undoable](img/sketchpad/storm-log-desktop.jpg)

### From words

Ask the map reads a storm without a model: "50 mm of rain over this area in 6 hours", "draw a
storm of 80 mm over 12 h", "rain 100 mm on Wales in 2 days with a runoff coefficient of 0.5",
"what if 120 mm of rain fell here". The action is `draw_storm` in `aquascope.map_commands`:
`rain_mm` (required), `duration_h` and `runoff_coefficient` (route_storm's own defaults, 6 h
and 0.3, when the words do not say), and the area: `area: "drawn"` ("this area", "here": the
area on the map, or the reader draws one), a `polygon`, a `bbox` or a `place` (its outline from
the gazetteer). A model may pick it; the package computes it.

![Ask the map: 80 mm of rain over this area in 12 hours](img/sketchpad/storm-ask-desktop.jpg)

On a phone the menu opens beside the tools, the card hides while the area is drawn, and after a
run it is a short sheet at the bottom over the time bar, with the area and the wave above it.

![Draw a storm on a phone: drawing the area](img/sketchpad/storm-draw-phone.jpg)
![Draw a storm on a phone: the area](img/sketchpad/storm-area-phone.jpg)
![Draw a storm on a phone: the wave and the card](img/sketchpad/storm-wave-phone.jpg)

### The check

The routing is Route it's, checked against GEOGLOWS's own hourly retrospective on four rivers
(NSE 0.89 to 0.99, the peak within 6 % and its timing within 5 h; the table is in
[Route a flood](explorer.md#the-check-against-geoglowss-own-simulation)). What the storm adds
is checked by its water:

- **The runoff is the rain times the coefficient times the area.** Over mid-Wales (the area in
  the pictures), 50 mm and 0.3 over 84 reaches whose unit catchments add up to 1,457.6 km² is
  0.3 × 0.050 m × 1,457.6 km² = 21.9 million m³, the runoff the package reports.
- **All of it leaves.** Routed for 185 h, the outlets pass 21.9 million m³: 100 % of the
  runoff (Muskingum keeps mass; `tests/test_routing.py` checks the same on a small network to
  1 %).

The runoff coefficient is not a model of the soil, and the routing is the flow the storm adds
to whatever the rivers already carry: modelled, for the shape and the timing of a wave, not a
forecast.

### The CLI and MCP

    aquascope route storm --polygon '[[-3.78,52.62],[-3.45,52.70],[-3.12,52.58],[-3.18,52.32],[-3.55,52.28],[-3.82,52.42]]' --rain 50 --duration 6 --runoff 0.3

prints the same message; `--json` adds every number, the hydrographs and the mass balance. The
MCP tool `route_storm` returns the same as JSON.

GEOGLOWS v2 routing configuration, CC BY 4.0; river network TDX-Hydro (NGA), CC BY-SA 4.0,
shown and not republished.
