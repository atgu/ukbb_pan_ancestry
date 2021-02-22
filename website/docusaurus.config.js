const math = require('remark-math')
const katex = require('rehype-katex')

module.exports = {
  title: 'Pan UKBB',
  tagline: 'Pan-ancestry genetic analysis of the UK Biobank',
  projectName: 'pan-ukbb',
  url: 'https://pan-dev.ukbb.broadinstitute.org',
  baseUrl: '/',
  favicon: 'img/favicon.ico',
  plugins: [
    [
      'remark-math',
      {
        id: 'remark-1',
      },
    ],
    [
      'rehype-katex',
      {
        id: 'katex-1',
      },
    ],
  ],
  stylesheets: ['https://cdn.jsdelivr.net/npm/katex@0.11.0/dist/katex.min.css'],
  themeConfig: {
    sidebarCollapsible: false,
    navbar: {
      // title: 'Pan UKBB',
      logo: {
        alt: 'My Site Logo',
        src: 'img/pan_ukbb_logo_v0.4-2.svg',
      },
      items: [
        {
          to: 'docs/summary',
          activeBasePath: 'docs',
          label: 'About',
          position: 'right',
        },
        { to: 'blog', label: 'News', position: 'right' },
        { to: 'downloads', label: 'Downloads', position: 'right' },
        { to: 'phenotypes', label: 'Phenotypes', position: 'right' },
        // {to: 'results', label: 'Results', position: 'right'},
        { to: 'team', label: 'Team', position: 'right' },
        { to: 'contact', label: 'Contact', position: 'right' },
        // {
        //   href: 'https://github.com/atgu/ukbb_pan_ancestry',
        //   label: 'GitHub',
        //   position: 'right',
        // },
      ],
    },
    footer: {
      style: 'dark',
      copyright: `Copyright © ${new Date().getFullYear()} Pan UKBB Team.`,
    },
    googleAnalytics: {
      trackingID: 'UA-169641392-1',
      anonymizeIP: true, // Should IPs be anonymized?
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          // It is recommended to set document id as docs home page (`docs/` path).
          sidebarPath: require.resolve('./sidebars.js'),
          remarkPlugins: [math],
          rehypePlugins: [katex],
        },
        blog: {
          showReadingTime: true,
          remarkPlugins: [math],
          rehypePlugins: [katex],
          // Please change this to your repo.
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
}
