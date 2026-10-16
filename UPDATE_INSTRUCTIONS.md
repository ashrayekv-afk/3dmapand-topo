# Update to NKPI 9.0.1

1. Use Windows **Extract All** on `NKPI_v901_Extract.zip`.
2. Copy/upload all extracted files and folders into the existing GitHub repository, replacing matching files. Keep the `components`, `nkpi`, `models`, and other folder structures intact.
3. Keep `kcn_dual_calculator_app.py` as the Streamlit main file. Reboot and refresh.
4. Verify `9.0.1-map-reader` in the footer.
5. Select Next patient / clear, then re-upload and process the original complete report ZIPs. Existing sessions retain the old reader results.
6. After calculation, open **Regional change & 3D corneal view → Current case → OD/OS → 3D corneal view**.

The new viewer works with one supported selected Belin visit, or with independently selectable earlier/latest visits. It needs readable elevation, BFS, and thickness data from the same eye and visit. It does not make a patient model from demo data or six scalar values. A failed optional spatial read does not change the calculator’s existing score behavior.

The release has 99 files and uses standard ZIP deflate compression. It replaces the obsolete standalone synthetic HTML preview with the patient 3D component. The in-app synthetic demo and numerical sample remain included. No source patient reports are included in the release.
