import { Datum } from "./types"

const capitalize = (text: string) => text.charAt(0).toUpperCase() + text.slice(1)

function prepareIcdPhenotype(phenotype: Datum) {
  let { description, description_more } = phenotype

  const icdDescriptionRe = /Date (\w+) first reported \((.+)\)/

  description = description || ''

  const match = description.match(icdDescriptionRe)

  if (match) {
    const icd10Code = match[1]
    let condition = match[2] || ''
    condition = capitalize(condition)
    description_more = `${description}; ${description_more}`
    description = `${icd10Code} ${condition}`
    return { ...phenotype, description, description_more }
  }

  return phenotype
}

function prepareCustomPhenotype(phenotype: Datum) {
  const { phenocode, modifier } = phenotype

  const customPhenoRe = /biogen/

  if (modifier) {
    return phenotype
  } else {
    const match = modifier.match(customPhenoRe)

    if (match) {
      return { ...phenotype, description: phenocode.replace(/_/g, ' ') }
    }

    return phenotype
  }

}

function prepareCodingDescriptionPhenotype(phenotype: Datum) {
  let { description, description_more } = phenotype
  const { coding_description } = phenotype

  const operationRe = /OPCS4|Non-cancer|Cancer code|Operation code|Treatment\/medication code/

  description = description || ''

  const match = description.match(operationRe)

  if (match) {
    description_more = `${description}; ${description_more}`
    description = capitalize(coding_description)
    return { ...phenotype, description, description_more }
  }

  return phenotype
}

const preparePrescriptionDescriptionPhenotype = (phenotype: Datum) => {
  const {trait_type, phenocode} = phenotype;
  if (trait_type === "prescriptions") {
    return {
      ...phenotype,
      description: phenocode,
    }
  } else {
    return phenotype
  }
}

const pipe = (...fns: ((phenotype: Datum) => Datum)[]) => (x: Datum) => fns.reduce((v, f) => f(v), x)

// Replace raw description field with custom values based on other fields.
// Adapted from UKBB browser.
export function processPhenotypeDescription(phenotype: Datum) {
  return pipe(
    prepareIcdPhenotype,
    prepareCodingDescriptionPhenotype,
    prepareCustomPhenotype,
    preparePrescriptionDescriptionPhenotype,
  )(phenotype)
}

