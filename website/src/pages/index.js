import React from 'react'
import classnames from 'classnames'
import Layout from '@theme/Layout'
import Link from '@docusaurus/Link'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import useBaseUrl from '@docusaurus/useBaseUrl'
import styles from './styles.module.css'

function Home() {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context
  return (
    <Layout
      title={`${siteConfig.title}`}
      description="Description will go into a meta tag in <head />"
    >
      <header className={classnames('hero', styles.heroBanner)}>
        <div className="container">
          <div className="text--center">
            <img className={styles.logo} src="img/pan_ukb_logo_v0.3-3.svg" alt={'pan ukbb title'} />
          </div>
          <p className="hero__subtitle">{siteConfig.tagline}</p>
        </div>
      </header>
      <main>
        <div className="container">
          <p class="hero__description">
            The goal of this project is to provide a resource to researchers that promotes more inclusive
            research practices, accelerates scientific discoveries, and improves the health of all people
            equitably. In genetics research, it is statistically necessary to study groups of individuals
            together with similar ancestries. In practice, this has meant that most previous research has
            excluded individuals with non-European ancestries. Here, we describe an effort to build a resource
            using one of the most widely accessed sources of genetic data, the UK Biobank, in a manner that is
            more inclusive than most previous efforts -- namely studying groups of individuals with diverse
            ancestries. The results of this research have a number of important limitations which should be
            carefully considered when researchers use this resource in their work and when they and others
            interpret subsequent findings.
          </p>

          <div className={styles.buttons}>
            <Link
              className={classnames('button button--outline button--secondary button--lg', styles.getStarted)}
              to={useBaseUrl('docs/')}
            >
              Learn more
            </Link>
          </div>
        </div>

      </main>
    </Layout>
  )
}

export default Home
