import React, { useState } from "react"
import { Accordion, AccordionDetails, AccordionSummary, FormControlLabel, FormGroup, FormLabel, makeStyles, Paper, Switch, Theme, Typography } from "@material-ui/core"
import { useCallback } from "react"
import { ColumnInstance, HeaderGroup, TableInstance, UseGlobalFiltersInstanceProps } from "react-table"
import { Datum } from "./types"
import { ExpandMore } from "@material-ui/icons"
import { ColumnGroupIndividualFilters } from "./ColumnGroupIndividualFilters"
import {  ColumnGroupFilterGroup, FilterDisplay, preventEventPropagation } from "./ColumnGroupFilterGroup";
import { GlobalFilter } from "./GlobalFilter"
import { commonPopulations, PopulationCode, populationColorMapping } from "./populations"
import { PopulationSwitch } from "./PopulationSwitch"

const useStyles = makeStyles((theme: Theme) => ({
  analysisAccordionTitle: {
    marginLeft: theme.spacing(6)
  },
  populationMetricsFilterGroup: {
    display: "grid",
    gridTemplateColumns: "repeat(3, 1fr)",
  }
}))

export type PopulationVisibility = {
  [K in PopulationCode]: boolean
}

export interface ColumnGroupVisibility {
  downloads: boolean,
  description: boolean,
  nCases: boolean,
  nControls: boolean,
  saigeHeritability: boolean,
  lambdaGc: boolean,
  md5: boolean,
}

enum AccordionName {
  Description,
  Analysis,
  Downloads,
  NCases,
  NControls,
  SaigeHeritability,
  LambdaGc,
  Md5,
}


interface Props {
  columnVisibilities: ColumnGroupVisibility
  setColumnVisibilities: React.Dispatch<React.SetStateAction<ColumnGroupVisibility>>
  populationMetricsVisibilities: PopulationVisibility
  setPopulationMetricsVisibilities: React.Dispatch<React.SetStateAction<PopulationVisibility>>
  columns: ColumnInstance<Datum>[]
  preGlobalFilteredRows: UseGlobalFiltersInstanceProps<Datum>["preGlobalFilteredFlatRows"]
  globalFilter: any
  setGlobalFilter: UseGlobalFiltersInstanceProps<Datum>["setGlobalFilter"]
}

export const PhenotypeFilters = (props: Props) => {
  const {
    columns, columnVisibilities, setColumnVisibilities,
    globalFilter, setGlobalFilter, preGlobalFilteredRows,
    populationMetricsVisibilities, setPopulationMetricsVisibilities,
  } = props;
  const classes = useStyles()

  const [expandedAccordion, setExpandedAccordion] = useState<AccordionName | undefined>(undefined)

  const handleColumnVisibilityChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    setColumnVisibilities({
      ...columnVisibilities, [event.target.name]: event.target.checked
  }) }, [columnVisibilities])

  const getAccordionChangeHandler = (accordionName: AccordionName) => (_: unknown, isExpanded: boolean) => setExpandedAccordion(isExpanded ? accordionName : undefined)

  const analysisFilter = (
    <Accordion expanded={expandedAccordion === AccordionName.Analysis} onChange={getAccordionChangeHandler(AccordionName.Analysis)}>
      <AccordionSummary expandIcon={<ExpandMore/>} className={classes.analysisAccordionTitle}>
        Analysis
      </AccordionSummary>
      <AccordionDetails>
        <FormGroup>
          <FormLabel>Filters</FormLabel>
          <ColumnGroupIndividualFilters
            columns={columns.find(col => col.columnGroupVisibilityAttribName === "analysis").columns}
          />
        </FormGroup>
      </AccordionDetails>
    </Accordion>
  )

  let descriptionFilterDisplay: FilterDisplay
  if (columnVisibilities.description) {
    descriptionFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={columns.find(col => col.columnGroupVisibilityAttribName === "description").columns}
        />
      )
    }
  } else {
    descriptionFilterDisplay = {showFilter: false}
  }

  const descriptionFilter = (
    <ColumnGroupFilterGroup
      visibilityControl={
        <Switch
          checked={columnVisibilities.description} name="description"
          onChange={handleColumnVisibilityChange}
        /> }
      label="Description"
      filterDisplay={descriptionFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.Description}
      onAccordionChange={getAccordionChangeHandler(AccordionName.Description)}
    />
  )

  const downloadsFilter = (
    <Paper>
      <AccordionSummary>
        <FormControlLabel
          onClick={preventEventPropagation}
          onFocus={preventEventPropagation}
          control={
            <Switch
              checked={columnVisibilities.downloads} name="downloads"
              onChange={handleColumnVisibilityChange}
            /> }
          label="Downloads"
        />
      </AccordionSummary>
    </Paper>
  )

  const handlePopulationMetricsVisibilityChange = (population: PopulationCode, isVisible: boolean) =>
    setPopulationMetricsVisibilities( currentState => ({...currentState, [population]: isVisible}))

  const populationMetricsVisibilitiesFilterElems = commonPopulations.map(pop => {
    const onChange = (event: React.ChangeEvent<HTMLInputElement>) =>
      handlePopulationMetricsVisibilityChange(pop, event.target.checked)
    const color = populationColorMapping.get(pop)
    return (
      <FormControlLabel key={pop}
        style={{color}}
        control={
          <PopulationSwitch
            population={pop}
            checked={populationMetricsVisibilities[pop]}
            onChange={onChange}
          />
        }
        label={pop}
      />
    )
  })

  const populationMetricsVibilitiesFilter = (
    <Accordion>
      <AccordionSummary expandIcon={<ExpandMore/>} className={classes.analysisAccordionTitle}>
        Per-Population Metrics
      </AccordionSummary>
      <AccordionDetails>
        <FormGroup className={classes.populationMetricsFilterGroup}>
          {populationMetricsVisibilitiesFilterElems}
        </FormGroup>
      </AccordionDetails>
    </Accordion>
  )


  let nCasesFilterDisplay: FilterDisplay
  if (columnVisibilities.nCases) {
    nCasesFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={columns.find(col => col.columnGroupVisibilityAttribName === "nCases").columns}
        />
      )
    }
  } else {
    nCasesFilterDisplay = {showFilter: false}
  }
  const nCasesFilter = (
    <ColumnGroupFilterGroup
      visibilityControl={
        <Switch
          checked={columnVisibilities.nCases} name="nCases"
          onChange={handleColumnVisibilityChange}
        /> }
      label="N Cases"
      filterDisplay={nCasesFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.NCases}
      onAccordionChange={getAccordionChangeHandler(AccordionName.NCases)}
    />
  )
  let nControlsFilterDisplay: FilterDisplay
  if (columnVisibilities.nControls) {
    nControlsFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={columns.find(col => col.columnGroupVisibilityAttribName === "nControls").columns}
        />
      )
    }
  } else {
    nControlsFilterDisplay = {showFilter: false}
  }
  const nControlsFilter = (
    <ColumnGroupFilterGroup
      visibilityControl={
        <Switch
          checked={columnVisibilities.nControls} name="nControls"
          onChange={handleColumnVisibilityChange}
        /> }
      label="N Controls"
      filterDisplay={nControlsFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.NControls}
      onAccordionChange={getAccordionChangeHandler(AccordionName.NControls)}
    />
  )
  let saigeHeritabilityFilterDisplay: FilterDisplay
  if (columnVisibilities.saigeHeritability) {
    saigeHeritabilityFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={[columns.find(col => col.columnGroupVisibilityAttribName === "saigeHeritability")]}
        />
      )
    }

  } else {
    saigeHeritabilityFilterDisplay = {showFilter: false}
  }
  const saigeHeritabilityFilter = (
    <ColumnGroupFilterGroup
      visibilityControl={
        <Switch
          checked={columnVisibilities.saigeHeritability} name="saigeHeritability"
          onChange={handleColumnVisibilityChange}
        /> }
      label="Saige Heritability"
      filterDisplay={saigeHeritabilityFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.SaigeHeritability}
      onAccordionChange={getAccordionChangeHandler(AccordionName.SaigeHeritability)}
    />
  )
  let lambdaGcFilterDisplay: FilterDisplay
  if (columnVisibilities.lambdaGc) {
    lambdaGcFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={columns.find(col => col.columnGroupVisibilityAttribName === "lambdaGc").columns}
        />
      )
    }

  } else {
    lambdaGcFilterDisplay = {showFilter: false}
  }
  const lambdaGcFilter = (
    <ColumnGroupFilterGroup
      visibilityControl={
        <Switch
          checked={columnVisibilities.lambdaGc} name="lambdaGc"
          onChange={handleColumnVisibilityChange}
        /> }
      label="Lambda GC"
      filterDisplay={lambdaGcFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.LambdaGc}
      onAccordionChange={getAccordionChangeHandler(AccordionName.LambdaGc)}
    />
  )
  const md5Filter = (
    <Paper>
      <AccordionSummary>
        <FormControlLabel
          onClick={preventEventPropagation}
          onFocus={preventEventPropagation}
          control={
            <Switch
              checked={columnVisibilities.md5} name="md5"
              onChange={handleColumnVisibilityChange}
            /> }
          label="MD5"
        />
      </AccordionSummary>
    </Paper>
  )

  const globalFilterElem = (
    <Paper>
      <AccordionDetails>
        <GlobalFilter
          preGlobalFilteredRows={preGlobalFilteredRows}
          globalFilter={globalFilter}
          setGlobalFillter={setGlobalFilter}
        />
      </AccordionDetails>
    </Paper>
  )


  return (
    <div>
      {globalFilterElem}
      {descriptionFilter}
      {analysisFilter}
      {downloadsFilter}
      {populationMetricsVibilitiesFilter}
      {nCasesFilter}
      {nControlsFilter}
      {saigeHeritabilityFilter}
      {lambdaGcFilter}
      {md5Filter}
    </div>
  )
}
