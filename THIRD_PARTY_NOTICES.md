# Third-Party Notices

## ChemDoodle Web Components

Files:

- `integrations/chemdoodle/ChemDoodleWeb.js`
- `integrations/chemdoodle/ChemDoodleWeb.css`

These files are ChemDoodle Web Components 11.0.0 from iChemLabs. Their headers
state that they are licensed under GNU GPL version 3 unless a separate
commercial iChemLabs license is used. Keep those headers intact when
redistributing the repository.

If you do not want GPLv3-covered files in a downstream distribution, remove the
two bundled ChemDoodle files and configure a local licensed copy with:

```bash
export CODEX_CHEM_CHEMDOODLE_WEB_DIR=/path/to/ChemDoodleWeb
```

## Ketcher npm packages

The Ketcher integration does not commit `node_modules/`. Its npm dependency
graph is resolved by `integrations/ketcher/package-lock.json` from these direct
dependencies:

- `@vitejs/plugin-react`
- `vite`
- `typescript`
- `react`
- `react-dom`
- `ketcher-react`
- `ketcher-standalone`

Review the package lock and upstream package metadata for exact transitive
versions and license terms.
