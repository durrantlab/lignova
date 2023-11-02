<h1 align="center">lignova</h1>

<h4 align="center">TODO</h4>

<h4 align="center" style="padding-bottom: 0.5em;"><a href="https://durrantlab.github.io/lignova">Documentation</a></h4>

<p align="center">
  <a href="https://github.com/psf/black" target="_blank">
    <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black style">
  </a>
  <a href="https://github.com/durrantlab/lignova/blob/main/.pre-commit-config.yaml" target="_blank">
    <img src="https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white" alt="License">
  </a>
  <a href="https://github.com/durrantlab/lignova/releases" target="_blank">
    <img src="https://img.shields.io/badge/%20%20%F0%9F%93%A6%F0%9F%9A%80-semantic--versions-e10079.svg" alt="License">
  </a>
  <a href="https://github.com/durrantlab/lignova/blob/main/LICENSE.md" target="_blank">
    <img src="https://img.shields.io/github/license/durrantlab/lignova" alt="License">
  </a>
</p>

TODO: Add in outline

## Deploying

A note to maintainers.

We use [bump-my-version](https://github.com/callowayproject/bump-my-version) to release a new version.
This will create a git tag that is used by [poetry-dynamic-version](https://github.com/mtkennerly/poetry-dynamic-versioning) to generate version strings and update `CHANGELOG.md`.

For example, to bump the `minor` version you would run the following command.

```bash
poetry run bump-my-version bump minor
```

After releasing a new version, you need to push and include all tags.

```bash
git push --follow-tags
```

## License

Code contained in this project is released under the MIT License as specified in [`LICENSE.md`][license].
This license grants you the freedom to use, modify, and distribute it as long as you include the original copyright notice contained in [`LICENSE.md`][license] and the following disclaimer.

> Portions of this code were incorporated and adapted with permission from [lignova](https://github.com/durrantlab/lignova) by durrantlab licensed under the [MIT License](https://github.com/durrantlab/lignova/blob/main/LICENSE.md).

[license]: https://github.com/durrantlab/lignova/blob/main/LICENSE.md
