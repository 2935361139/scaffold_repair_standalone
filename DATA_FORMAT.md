# Optional trajectory contract

Endpoint-only use needs no trajectory. Trajectories enable observed timing, not proof of an upstream cause.

For a generic manifest sample, add `trajectory_npz` and `trajectory_metadata`.
All trajectory atoms use the generated SDF's atom order, including fixed atoms.
Coordinates are already in the protein/world frame, in angstroms. No centering or alignment is performed.

Metadata example:

```json
{
  "coordinate_frame": "protein_world_angstrom",
  "atom_vocabulary": ["C", "N", "O", "S"],
  "terminal_stream": "sampler",
  "streams": {
    "sampler": {
      "indices": [0, 1, 2],
      "times": [0.0, 0.5, 1.0],
      "stages": ["state", "state", "final"]
    }
  }
}
```

NPZ arrays for each stream named `sampler`:

- `sampler_coords`: finite float array `(T,N,3)`.
- `sampler_atom_classes`: integer array `(T,N)` indexing `atom_vocabulary`.
- `sampler_valid`: boolean array `(T,N)`.

Frame indices must be consecutive from zero. Native time can be null where unavailable.
The terminal stream's last coordinates must match the high precision endpoint (absolute tolerance `2e-5 Å`), with the same elements and all retained atoms active.
During intermediate frames a missing/different element is recorded as unknown for that final clash pair.
The numerical ambiguity band `1e-5 Å` in trajectory observations is unchanged from the previous audit; it does not change the endpoint PoseBusters threshold.

Variable atom counts, atom permutations, and different coordinate frames must be normalized by the model adapter before using this schema. They are not guessed from the SDF.

FLOWR records are imported by their own adapter. It validates the saved molecule order, fixed mask, reference, protein hash and original endpoint labels, then converts the recorded trajectory to the protein frame.
Where original `inputs.pt` is unavailable, the adapter recovers the translation only if the float32 rounding intervals identify one unique translation reproducing every terminal coordinate exactly. Otherwise it stops. It does not fit an approximate alignment.

## Output evidence

`generated/samples/<id>/coordinates.npz` retains all normalized original coordinates and conditioning targets. `raw.sdf`, `reference.sdf`, `protein.pdb`, and optional `trajectory.npz` make the checked dataset self-contained.

`generated/integrity.json` checks this one batch. It is not a history log and does not exclude reused ligands. This protects separation of checking and intervention: a modified molecule, mapping, protein or method cannot silently enter the old experiment.

FLOWR raw-graph fixed-identity labels for chemically invalid/undecodable outputs are preserved from the checksum-verified previous audit and explicitly labeled `previous_checksum_verified_raw_graph_audit`; they are not claimed to be newly inferred from an empty SDF. Such samples are unsupported by the coordinate repair bank. For readable valid molecules, the fixed-identity check is recomputed from the mapped graph.

External SDF stereochemical annotations are preserved. In the FLOWR record adapter only, reader-inferred 3D chirality is removed to match the original tensor decoder, which did not assign those tags; this is necessary for exact historical output comparison.

Classification never drops a selected SDF record silently. Unreadable/invalid generated molecular records are retained as failed/unsupported. Missing files, incorrect mappings, evaluator exceptions, changed evidence, or incomplete baselines stop a stage with an error; they are not counted as successful molecular evaluations.
