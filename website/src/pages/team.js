import React from 'react'
import classnames from 'classnames'
import Layout from '@theme/Layout'
import Link from '@docusaurus/Link'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import useBaseUrl from '@docusaurus/useBaseUrl'
import styles from './styles.module.css'

function Contributors() {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context
  return (
    <Layout title={`${siteConfig.title}`} description="Team">
      <header className={classnames('hero', styles.heroBanner)}>
        <div className="container">
          <h1>Team</h1>
        </div>
      </header>
      <main>
        <div className="container">
          <p>
            <b>Analysis</b>
            <ul>
              <li>Konrad Karczewski</li>
              <li>Elizabeth Atkinson</li>
              <li>Masahiro Kanai</li>
              <li>Nikolas Baya</li>
              <li>Patrick Turley</li>
              <li>Gopal Sarma</li>
              <li>Raymond Walters</li>
              <li>Duncan Palmer</li>
              <li>Nathan Cheng</li>
              <li>Rahul Gupta</li>
              <li>Wei Zhou</li>
              <li>Alicia Martin</li>
            </ul>
            <b>Ethics and communication</b>
            <ul>
              <li>Shawneequa Callier</li>
              <li>Patrick Turley</li>
              <li>Elizabeth Atkinson</li>
              <li>Alicia Martin</li>
            </ul>
            <b>Website</b>
            <ul>
              <li>Matthew Solomonson</li>
              <li>Konrad Karczewski</li>
            </ul>
            <b>Data coordination</b>
            <ul>
              <li>Sam Bryant</li>
              <li>Claire Churchhouse</li>
              <li>Caroline Cusick</li>
            </ul>
            <b>Hail</b>
            <ul>
              <li>Jacqueline Goldstein</li>
              <li>Daniel King</li>
              <li>Cotton Seed</li>
            </ul>
            <b>Coordination and Leadership</b>
            <ul>
              <li>Alicia Martin</li>
              <li>Hilary Finucane</li>
              <li>Benjamin Neale</li>
              <li>Mark Daly</li>
            </ul>
          </p>
        </div>
      </main>
    </Layout>
  )
}

export default Contributors
