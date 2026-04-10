# Callable Domains

The advanced surfaces in `paired_sc` are API domains, not sidecar
features hidden behind packaging switches.

Supported domains:

- `liana`
- `magic`
- `trajectory`
- `latent`
- `regulatory`

Each domain exposes:

- a Python callable
- a CLI subcommand
- structured status reporting
- isolated outputs under `results/domains/<domain>/`

If a domain dependency or required input is unavailable, the domain returns a
structured skipped status instead of crashing the whole run.

