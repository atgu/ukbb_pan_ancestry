# Website

This website is built using [Docusaurus 2](https://v2.docusaurus.io/), a modern static website generator.

### Installation

```
$ yarn
```

### Local Development

```
$ yarn start
```

This command starts a local development server and open up a browser window. Most changes are reflected live without having to restart the server.

### Build

```
$ yarn build
```

This command generates static content into the `build` directory and can be served using any static contents hosting service.

### Deployment

Fix a bug in `@docusaurus/core@^2.0.0-alpha.56`.

`npm run patch-route-bug`

Build and upload build to storage bucket.

`npm run deploy-dev`

### Manifest table data update

The manifest table shown in the UI is read from `website/src/data.json`.
To update this file, go to the online Google Sheet and download is as CSV (File > Download > CSV) and convert the download CSV to JSON, for example, using [this online tool](https://csvjson.com/csv2json).
