# The sketchpad: draw on the map, and a small model recomputes

The sketchpad is part of [make the map compute](https://github.com/Rekin226/aquascope/issues/664): you draw on the
globe, a small, checked model runs in your browser, and the answer is drawn back on the map. Each tool is one plain
function in `aquascope.sketch` (experimental), with a CLI verb, an MCP tool and an Ask the map action over the same
function. Every drawing goes in the action log on the map, so one click undoes it.

## Drop a well

"What happens to the water table if I pump here?" The well tool in the map tools (the pump mark) drops a pumping well
where you click. Its drawdown cone is drawn round it: contour lines at 5 cm (dashed, the cone's edge), 10 cm, 20 cm,
0.5 m, 1, 2, 5 m and deeper, over a fill that darkens with depth. Drop more wells and their cones add. Drag a well and
the cone follows as you drag; change a well's rate, the aquifer's transmissivity and storativity, or the days of
pumping, and it is computed again. **Play the cone growing** takes the time bar and plays the cone from 6 hours to the
days you chose.

A well's card gives its drawdown at its face (at the well's radius, 0.15 m by default, with no well loss), its own
share and what the other wells add there, and how far its 5 cm line reaches. A contour's card gives the drawdown
where you clicked and the nearest well. On a phone the panel is a sheet at the bottom, and a card takes its place
while it is open.

![Two wells near Lyon, after 30 days: the cones add](img/sketchpad/well-two-desktop.jpg)

![The cone growing on the time bar, after 3 days](img/sketchpad/well-growing-desktop.jpg)

**With your MODFLOW model.** When MODFLOW heads are on the map ([MODFLOW heads on the map](modflow.md)) and a well is
dropped inside the model, its card names the model's cell under it, the cell's size and the layer shown, and what the
model's well cell would read instead of the well's own drawdown: the Theis drawdown at Peaceman's equivalent radius,
0.2 times the cell's side (Peaceman, 1978). A coarse cell averages the steep part of the cone; the card says by how
much.

![A well over the Freyberg sample: the model's cell and Peaceman's radius](img/sketchpad/well-modflow-desktop.jpg)

### The model

Theis (1935): the drawdown at a distance `r` from a well pumping `Q` for a time `t` in a confined aquifer of
transmissivity `T` and storativity `S` is

```
s = Q / (4 π T) · W(u),    u = r² S / (4 T t),    W(u) = E1(u)
```

Several wells add by superposition. The well function is `aquascope.groundwater.aquifer`'s, the engine behind
`aquascope groundwater --analysis theis` and the `aquifer_drawdown` tool, so there is one Theis in the package. The
cone is traced by marching squares on a grid that is even across the cone and closes in on every well geometrically,
so the tight circles near a well are as true as the wide ones far out. Distances are measured on a local plane around
the wells. The cone is cut 100 km from the wells.

**What it assumes**, and the panel says so: a confined aquifer, homogeneous, isotropic and of infinite extent; wells
that fully penetrate it and pump at a constant rate; water released from storage at once. No recharge, no rivers or
boundaries, no leakage, no delayed yield of a water table. A screening answer, not a well-field design.

**The defaults** are a municipal well in a confined sand aquifer: 1,000 m³/day (11.6 L/s), T = 500 m²/day,
S = 1e-4, read after 30 days. The typical ranges are from Heath (1983), USGS Water-Supply Paper 2220: the
transmissivity of principal aquifers ranges from less than 1 to 100,000 m²/day (p. 26), and the storage coefficient of
most confined aquifers from about 1e-5 to 1e-3 (p. 28). The panel's sliders reach past them, and say so when you do.

### The check

`aquascope well --check` recomputes these, and `tests/test_sketch_drawdown.py` holds them:

| Reference | Printed | Computed |
| --- | --- | --- |
| Fetter (2001), *Applied Hydrogeology*, 4th ed., the worked example of the Theis equation: K = 14.9 m/day over 20.1 m (T = 299 m²/day as printed), S = 0.0051, Q = 2725 m³/day, r = 7.0 m, t = 1 day | u = 0.00021, W(u) = 7.94, drawdown 5.7 m | u = 0.000209, drawdown 5.73 m |
| Lohman (1972), *Ground-water hydraulics*, USGS Professional Paper 708, Table 5, p. 16: W(u) for u = N × 10ᵏ, rows N = 1.0 and N = 2.0, k from -15 to 0 (32 values) | four decimals | every value within 1e-4 |
| Superposition: two equal wells 1 km apart | equal at mirror points, the same both ways | 2.4021 / 2.4021 m; 1.1290 / 1.1290 m |

Fetter reads W(u) from the table at u = 2 × 10⁻⁴ (7.94) where the inputs give u = 2.09 × 10⁻⁴ (W = 7.90); the book's
5.7 m and the code's 5.73 m agree to the two figures printed. PP 708 is at
<https://pubs.usgs.gov/pp/0708/report.pdf>; the superposition check also confirms that the wells pumped together give
the wells pumped alone, added.

### From words, a CLI and an agent

Ask the map reads "drop a well here", "drop a well at Lyon pumping 2000 m3/day" and "drill a borehole near Lille
pumping 20 l/s for 90 days" into the `drop_well` action without a model (a rate in L/s is turned into m³/day); the
place is looked up and the map computes the cone. The same function from a terminal and for an agent:

```bash
aquascope well --at 4.88 45.76 --at 4.92 45.765 --rate 1000 --rate 3000 --days 90 --geojson cone.geojson
aquascope well --check
```

![Ask the map: a borehole near Vaulx-en-Velin, 20 L/s for 90 days, in the undo log](img/sketchpad/well-ask-desktop.jpg)

The MCP tool `drop_well` takes one well (`lon`, `lat`, `rate_m3_day`) or `wells`, and returns the same numbers and the
contours in lon/lat.

![Drop a well on a phone](img/sketchpad/well-phone.jpg)

### References

- Theis, C.V. (1935). The relation between the lowering of the piezometric surface and the rate and duration of
  discharge of a well using ground-water storage. *Transactions, American Geophysical Union*, 16, 519-524.
- Fetter, C.W. (2001). *Applied Hydrogeology*, 4th edition. Prentice Hall.
- Lohman, S.W. (1972). *Ground-water hydraulics*. USGS Professional Paper 708.
- Heath, R.C. (1983). *Basic ground-water hydrology*. USGS Water-Supply Paper 2220.
- Peaceman, D.W. (1978). Interpretation of well-block pressures in numerical reservoir simulation. *SPE Journal*,
  18(3), 183-194. doi:10.2118/6893-PA.
