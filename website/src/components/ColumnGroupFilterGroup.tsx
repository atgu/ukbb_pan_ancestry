import React from "react"
import { Accordion, AccordionDetails, AccordionProps, AccordionSummary, FormControlLabel, FormGroup, FormLabel } from "@material-ui/core"
import { ExpandMore } from "@material-ui/icons"

export const preventEventPropagation = (event: React.MouseEvent | React.FocusEvent) => event.stopPropagation()

export type FilterDisplay = {
  showFilter: false
} | {
  showFilter: true
  filterElems: React.ReactNode
}

interface Props {
  visibilityControl: React.ReactNode
  label: string
  filterDisplay: FilterDisplay
  isAccordionExpanded: AccordionProps["expanded"]
  onAccordionChange: AccordionProps["onChange"]
}

export const ColumnGroupFilterGroup = (props: Props) => {
  const {visibilityControl, filterDisplay, label, isAccordionExpanded, onAccordionChange} = props
  let accordionDetails: React.ReactNode
  if (filterDisplay.showFilter) {
    accordionDetails = (
      <FormGroup>
        <FormLabel>Filters</FormLabel>
        {filterDisplay.filterElems}
      </FormGroup>
    )
  } else {
    accordionDetails = (
      <FormGroup>
        <FormLabel>Filters only available when columns are shown</FormLabel>
      </FormGroup>
    )
  }
  return (
    <Accordion expanded={isAccordionExpanded} onChange={onAccordionChange}>
      <AccordionSummary expandIcon={<ExpandMore/>} >
        <FormControlLabel
          onClick={preventEventPropagation}
          onFocus={preventEventPropagation}
          control={visibilityControl}
          label={label}
        />
      </AccordionSummary>
      <AccordionDetails>
        {accordionDetails}
      </AccordionDetails>
    </Accordion>
  )
}
