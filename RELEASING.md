# Releasing Twind to PyPI

Releases are published automatically by the GitHub Actions workflow at
`.github/workflows/publish.yml`. The workflow triggers whenever a tag matching
`v*` is pushed to the repository (it can also be run manually from the Actions
tab).

## One-time setup

The workflow authenticates to PyPI with an API token stored as the GitHub
secret `PYPI_API_TOKEN`. To rotate or recreate it:

1. Create a project-scoped token at
   https://pypi.org/manage/account/token/ (scope: project `twind`).
2. In the GitHub repo: `Settings` → `Secrets and variables` → `Actions` →
   update or add `PYPI_API_TOKEN`.

## Cutting a new release

1. Make and commit your code changes on `master`.
2. Bump the version in `twind/__version__.py` (follow semver: patch for bug
   fixes, minor for backwards-compatible features, major for breaking changes).
3. Commit the bump:
   ```bash
   git add twind/__version__.py
   git commit -m "Bump version to X.Y.Z"
   ```
4. Tag the commit and push both the branch and the tag:
   ```bash
   git tag vX.Y.Z
   git push origin master
   git push origin vX.Y.Z
   ```
5. Watch the run at https://github.com/changgoo/Twind/actions. On success, the
   new version appears at https://pypi.org/project/twind/.

## Notes

- PyPI rejects duplicate versions, so forgetting to bump
  `twind/__version__.py` will cause the publish step to fail rather than
  overwrite the existing release.
- Tag names must start with `v` (e.g. `v1.1.3`) to match the workflow trigger.
- To re-run a failed publish without changing code, use `Run workflow` on the
  Actions page (`workflow_dispatch`).
