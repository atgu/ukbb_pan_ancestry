module.exports = {
  title: 'Pan UKBB',
  tagline: 'Pan-ancestry genetic analysis of the UK Biobank',
  url: 'https://pan-dev.ukbb.broadinstitute.org',
  baseUrl: '/',
  favicon: 'img/favicon.ico',
  organizationName: 'atgu', // Usually your GitHub org/user name.
  projectName: 'ukbb_pan_ancestry', // Usually your repo name.
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
          label: 'FAQ',
          position: 'left',
        },
        {
          to: 'docs/downloads',
          activeBasePath: 'docs',
          label: 'Downloads',
          position: 'left',
        },
        {to: 'blog', label: 'News', position: 'left'},
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
          // Please change this to your repo.
          editUrl:
            'https://github.com/atgu/ukbb_pan_ancestry/edit/master/website',
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/atgu/ukbb_pan_ancestry/edit/master/website/blog',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
