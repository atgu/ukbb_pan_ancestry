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
import { ColumnGroupVisibility, PerPopulationMetricsVisibility, ColumnGroupName } from "./phenotypesReducer"

const useStyles = makeStyles((theme: Theme) => ({
  analysisAccordionTitle: {
    marginLeft: theme.spacing(6)
  },
  populationMetricsFilterGroup: {
    display: "grid",
    gridTemplateColumns: "repeat(3, 1fr)",
  }
}))


enum AccordionName {
  Description,
  Analysis,
  Downloads,
  NCases,
  NControls,
  SaigeHeritability,
  LambdaGc,
  Md5,
  PerPopulationMetrics
}


interface Props {
  columnVisibilities: ColumnGroupVisibility
  setColumnVisibilities: (population: ColumnGroupName, isVisible: boolean) => void
  populationMetricsVisibilities: PerPopulationMetricsVisibility
  setPopulationMetricsVisibilities: (population: PopulationCode, isVisible: boolean) => void
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
    setColumnVisibilities(event.target.name as ColumnGroupName, event.target.checked)
  }, [columnVisibilities])

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
            columns={columns.find(col => col.columnGroupName === ColumnGroupName.Analysis).columns}
          />
        </FormGroup>
      </AccordionDetails>
    </Accordion>
  )

  let descriptionFilterDisplay: FilterDisplay
  if (columnVisibilities[ColumnGroupName.Description]) {
    descriptionFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={columns.find(col => col.columnGroupName === ColumnGroupName.Description).columns}
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
          checked={columnVisibilities[ColumnGroupName.Description]} name={ColumnGroupName.Description}
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
              checked={columnVisibilities[ColumnGroupName.Downloads]} name={ColumnGroupName.Downloads}
              onChange={handleColumnVisibilityChange}
            /> }
          label="Downloads"
        />
      </AccordionSummary>
    </Paper>
  )


  const populationMetricsVisibilitiesFilterElems = commonPopulations.map(pop => {
    const onChange = (event: React.ChangeEvent<HTMLInputElement>) =>
      setPopulationMetricsVisibilities(pop, event.target.checked)
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
    <Accordion expanded={expandedAccordion === AccordionName.PerPopulationMetrics} onChange={getAccordionChangeHandler(AccordionName.PerPopulationMetrics)}>
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
  if (columnVisibilities[ColumnGroupName.NCases]) {
    nCasesFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={columns.find(col => col.columnGroupName === ColumnGroupName.NCases).columns}
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
          checked={columnVisibilities[ColumnGroupName.NCases]} name={ColumnGroupName.NCases}
          onChange={handleColumnVisibilityChange}
        /> }
      label="N Cases"
      filterDisplay={nCasesFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.NCases}
      onAccordionChange={getAccordionChangeHandler(AccordionName.NCases)}
    />
  )
  let nControlsFilterDisplay: FilterDisplay
  if (columnVisibilities[ColumnGroupName.NControls]) {
    nControlsFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={[columns.find(col => col.columnGroupName === ColumnGroupName.NControls)]}
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
          checked={columnVisibilities[ColumnGroupName.NControls]} name={ColumnGroupName.NControls}
          onChange={handleColumnVisibilityChange}
        /> }
      label="N Controls"
      filterDisplay={nControlsFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.NControls}
      onAccordionChange={getAccordionChangeHandler(AccordionName.NControls)}
    />
  )
  let saigeHeritabilityFilterDisplay: FilterDisplay
  if (columnVisibilities[ColumnGroupName.SaigeHeritability]) {
    saigeHeritabilityFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={[columns.find(col => col.columnGroupName === ColumnGroupName.SaigeHeritability)]}
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
          checked={columnVisibilities[ColumnGroupName.SaigeHeritability]} name={ColumnGroupName.SaigeHeritability}
          onChange={handleColumnVisibilityChange}
        /> }
      label="Saige Heritability"
      filterDisplay={saigeHeritabilityFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.SaigeHeritability}
      onAccordionChange={getAccordionChangeHandler(AccordionName.SaigeHeritability)}
    />
  )
  let lambdaGcFilterDisplay: FilterDisplay
  if (columnVisibilities[ColumnGroupName.LambdaGc]) {
    lambdaGcFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={columns.find(col => col.columnGroupName === ColumnGroupName.LambdaGc).columns}
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
          checked={columnVisibilities[ColumnGroupName.LambdaGc]} name={ColumnGroupName.LambdaGc}
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
              checked={columnVisibilities[ColumnGroupName.Md5]} name={ColumnGroupName.Md5}
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
