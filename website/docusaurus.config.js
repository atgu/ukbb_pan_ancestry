module.exports = {
  title: 'Pan UKBB',
  tagline: 'Pan-ancestry genetic analysis of the UK Biobank',
  url: 'https://pan-dev.ukbb.broadinstitute.org',
  baseUrl: '/',
  favicon: 'img/favicon.ico',
  themeConfig: {
    sidebarCollapsible: false,
    navbar: {
      title: 'Pan UKBB',
      // logo: {
      //   alt: 'My Site Logo',
      //   src: 'img/logo.svg',
      // },
      links: [
        {
          to: 'docs/background',
          activeBasePath: 'docs',
          label: 'About',
          position: 'left',
        },
        {to: 'blog', label: 'News', position: 'left'},
        {to: 'downloads', label: 'Downloads', position: 'left'},
        {to: 'results', label: 'Results', position: 'left'},
        {
          href: 'https://github.com/atgu/ukbb_pan_ancestry',
          label: 'GitHub',
          position: 'right',
        },
      ],
    },
    footer: {
      style: 'dark',
      copyright: `Copyright © ${new Date().getFullYear()} Pan UKBB Team.`,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          // It is recommended to set document id as docs home page (`docs/` path).
          homePageId: 'background',
          sidebarPath: require.resolve('./sidebars.js'),
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
