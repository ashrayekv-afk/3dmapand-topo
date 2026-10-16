# Patient corneal shape — version 9.0.1

## Two separate sources

**Report-image estimate:** the standard anterior elevation map and printed floating BFS radius from the selected Belin–Ambrósio page, together with that page's thickness map. Front shape is a centered, untilted sphere plus signed elevation, converted from µm to mm. The printed floating BFS translation and tilt are unavailable; this reconstruction explicitly assumes their values. Enhanced/exclusion and difference maps are not substituted. Curvature colors and NKPI scalars never generate height.

A local weighted plane (0.6-mm radius) smooths the elevation and thickness. At least six valid neighboring samples and 50% local readable support are required, with a maximum nearest-sample distance of two sample steps (or 0.25 mm for very fine grids). Large missing areas remain open. A model requires at least 50% coverage of the sampled central region. The inner surface is an offset along the reconstructed anterior normal by image-derived pachymetry. This is a modeling approximation, not recovery of independently measured posterior geometry. The displayed range is the modeled thickness range, not the device's minimum thickness.

**Native height reconstruction:** matching `kind: measured` / `source: device_numeric_export` anterior and optional posterior `z_mm` grids pass the existing strict `CornealData` schema. Declared physical axes, laterality and z direction are respected. A supplied device thickness map colors the mesh; the renderer does not infer thickness from raw vertical front/back separation. Missing cells and reversed front/back ordering are omitted. Synthetic fixtures cannot enter the patient route.

## Association and display

The selected eye, visit and input snapshot must match the imported source binding. Input changes, source changes, removal of an eye and next-patient reset sever the link. The view never silently uses the other eye or a simulated cone. Earlier/latest shapes are independent views, not a registered height subtraction. OD and OS keep the native report orientation (OD temporal at screen left, OS temporal at screen right in front view).

All physical axes use the same scale. There is no vertical exaggeration. The native or estimated source label is embedded in exported PNGs. The component uses local Canvas rendering and the existing local OCR models; there is no new API or remote AI reader. No report images, calibration grids or geometry are sent to the score server.

## Native export contract

The existing `nkpi-corneal-surface-1` schema requires:
- eye and visit (`OD`/`OS`, `first`/`latest`), measured source and units;
- finite, increasing x/y axes in mm and explicit orientation;
- anterior `z_mm`, optional posterior `z_mm`, optional device `thickness_um`;
- null for missing samples; enough complete cells for rendering.

Use the existing `nkpi-surface-index.json` association workflow in `docs/AUTO_DICOM_IMPORT.md` and numerical conversion documentation. Image-derived estimates have their own internal `nkpi-image-shape-1` contract and are never labeled as device exports.

## Basis and limitations

Reference-sphere radius and placement affect interpretation of elevation maps: https://pmc.ncbi.nlm.nih.gov/articles/PMC4142625/ . The report does not supply all of that pose information. Therefore this supplement is an approximate visualization, not an exact anatomical replica, a validated diagnostic surface, a surgical planning model, or evidence that a region will progress. The calculator's scalar inputs and NKPI are unchanged.

## Diameter confirmation and sample coverage

Version 9.0.1 uses expanded tag crops, cross-crop OCR, numeral-shape recognition and an isolated-digit retry. A valid printed plot size and unit are required. A BFS fitting diameter is never used as the plot diameter. The shape view uses the updated image sampler, including explicitly tracked local-color recovery and rejection of isolated colored-marker artifacts. Modeled coverage includes the further local smoothing/interpolation described above; it is not a direct-data or accuracy percentage.
