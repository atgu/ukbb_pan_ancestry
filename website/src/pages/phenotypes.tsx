import React from 'react'
import Layout from '@theme/Layout'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import {  PhenotypesPageContentWrapper } from "../components/PhenotypesPageContentWrapper";
import "core-js/stable";
import "regenerator-runtime/runtime";

const Phenotypes = () => {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context
  return (
    <Layout title={`${siteConfig.title}`} description="Phenotypes">
      <PhenotypesPageContentWrapper/>
    </Layout>
  )
}

export default Phenotypes
