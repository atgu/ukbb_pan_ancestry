import React from 'react'
import Layout from '@theme/Layout'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import {  PhenotypesPageContentWrapper } from "../components/PhenotypesPageContentWrapper";
import {  docusaurusLayoutWrapperClassName } from "../components/PhenotypesPageContent";

import "regenerator-runtime/runtime.js";
import "core-js/stable";

const Phenotypes = () => {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context
  return (
    <Layout
      title={`${siteConfig.title}`}
      description="Phenotypes"
      wrapperClassName={docusaurusLayoutWrapperClassName}
    >
      <PhenotypesPageContentWrapper/>
    </Layout>
  )
}

export default Phenotypes
