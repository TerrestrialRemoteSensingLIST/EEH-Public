"""Supporting tools for the EEH2 pipeline.

- ``download_sample.py`` — fetch ECOSTRESS and all ancillary inputs
- ``h5_to_cog.py`` — convert algorithm HDF5 outputs to Cloud Optimized GeoTIFF

Both are runnable as scripts; ``h5_to_cog`` is also imported by
``run_pipeline.py`` for the ``--cog`` step.
"""
