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
    <Layout title={`${siteConfig.title}`} description="Description will go into a meta tag in <head />">
      <header className={classnames('hero', styles.heroBanner)}>
        <div className="container">
          <div className="text--center">
            <img
              className={styles.logo}
              id="pan_ukbb_logo"
              src="img/pan_ukbb_logo_v0.4-2.svg"
              alt={'pan ukbb title'}
            />
          </div>
          <p className="hero__subtitle">{siteConfig.tagline}</p>
        </div>
      </header>
      <main>
        <div className="container">
          <p className="hero__description">
            The UK Biobank is a collection of a half million individuals with paired genetic and phenotype
            information that has been enormously valuable in studies of genetic etiology for common diseases
            and traits. However, most genome-wide analyses of this dataset use only the European ancestry
            individuals. Analyzing a more inclusive and diverse dataset increases power and improves the
            potential for discovery. Here, we present a multi-ancestry analysis of 7,221 phenotypes, across 6
            continental ancestry groups, for a total of 16,119 genome-wide association studies. We release
            these summary statistics freely to the community ahead of publication.
          </p>
          <p className="hero__description">
            We acknowledge that genetic studies in diverse-ancestry individuals are culturally sensitive. We
            therefore regularly engaged with and sought feedback from a local affinity group for the
            advancement of ethnic minorities, several scholars and clinicians with experience working with
            minority groups, and a bioethics expert. They helped direct and contextualize this research,
            discussed actions to maximize benefits and minimize risks, and provided feedback on{' '}
            <Link to={useBaseUrl('docs/summary')}>an FAQ page</Link> describing our study, its risks,
            benefits, and limitations. We intend for this resource to promote more inclusive research
            practices, accelerate novel scientific discoveries that are more generalizable, and improve the
            health of all people equitably.
          </p>
          <p className="hero__description">
            Together, these efforts vastly extend the available resources for interpretation of
            disease-associated variants across diverse populations and highlight the importance of increasing
            diversity in genetic studies.
          </p>

          <div className={styles.buttons}>
            <Link
              className={classnames('button button--outline button--secondary button--lg', styles.getStarted)}
              to={useBaseUrl('docs/background')}
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
