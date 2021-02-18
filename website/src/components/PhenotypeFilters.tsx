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
import { ColumnGroupVisibility, PerPopulationMetricsVisibility, ColumnGroupName, PerPopulationRangeFilter, RangeFilterMetric, RangeFilterValue } from "./phenotypesReducer"
import {  NumberRangeColumnFilter } from "./NewNumberRangeColumnFilter";


const useStyles = makeStyles((theme: Theme) => ({
  analysisAccordionTitle: {
    marginLeft: theme.spacing(6)
  },
  populationMetricsFilterGroup: {
    display: "grid",
    gridTemplateColumns: "repeat(3, 1fr)",
  }
}))

export type PerPopulationExtremums = {
  [K in PopulationCode]: {min: number, max: number}
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
  nCasesFilters: PerPopulationRangeFilter
  nCasesPerPopulationExtremums: PerPopulationExtremums
  nControlsFilters: PerPopulationRangeFilter
  nControlsPerPopulationExtremums: PerPopulationExtremums
  saigeHeritabilityFilters: PerPopulationRangeFilter
  saigeHeritabilityPerPopulationExtremums: PerPopulationExtremums
  lambdaGcFilters: PerPopulationRangeFilter
  lambdaGcPerPopulationExtremums: PerPopulationExtremums
  disableFilterOnePopulation: (args: {metric: RangeFilterMetric, population: PopulationCode}) => void
  updateFilterOnePopulation: (args: {metric: RangeFilterMetric, population: PopulationCode, min: number, max: number}) => void
}

export const PhenotypeFilters = (props: Props) => {
  const {
    columns, columnVisibilities, setColumnVisibilities,
    globalFilter, setGlobalFilter, preGlobalFilteredRows,
    populationMetricsVisibilities, setPopulationMetricsVisibilities,
    updateFilterOnePopulation, disableFilterOnePopulation,
    nCasesFilters, nCasesPerPopulationExtremums,
    nControlsFilters, nControlsPerPopulationExtremums,
    saigeHeritabilityFilters, saigeHeritabilityPerPopulationExtremums,
    lambdaGcFilters, lambdaGcPerPopulationExtremums,
  } = props;
  const classes = useStyles()

  const [expandedAccordion, setExpandedAccordion] = useState<AccordionName | undefined>(undefined)

  const handleColumnVisibilityChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    setColumnVisibilities(event.target.name as ColumnGroupName, event.target.checked)
  }, [columnVisibilities])

  const getAccordionChangeHandler = (accordionName: AccordionName) => (_: unknown, isExpanded: boolean) => setExpandedAccordion(isExpanded ? accordionName : undefined)

  let analysisFilterDisplay: FilterDisplay
  if (columnVisibilities[ColumnGroupName.Analysis]) {
    analysisFilterDisplay = {
      showFilter: true,
      filters: (
        <ColumnGroupIndividualFilters
          columns={columns.find(col => col.columnGroupName === ColumnGroupName.Analysis).columns}
        />
      )
    }
  } else {
    analysisFilterDisplay = {showFilter: false}
  }
  const analysisFilter = (
    <ColumnGroupFilterGroup
      visibilityControl={
        <Switch
          checked={columnVisibilities[ColumnGroupName.Analysis]} name={ColumnGroupName.Analysis}
          onChange={handleColumnVisibilityChange}
        /> }
      label="Analysis"
      filterDisplay={analysisFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.Analysis}
      onAccordionChange={getAccordionChangeHandler(AccordionName.Analysis)}
    />

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
    const nCasesPopulationFilterElems = Object.entries(nCasesFilters).map(([population, filterValue]) => {
      return (
        <NumberRangeColumnFilter key={population}
          label={population}
          population={population as PopulationCode}
          globalMinValue={nCasesPerPopulationExtremums[population as PopulationCode].min}
          globalMaxValue={nCasesPerPopulationExtremums[population as PopulationCode].max}
          metric={RangeFilterMetric.NCases}
          filterValue={filterValue}
          disableFilter={disableFilterOnePopulation}
          updateFilter={updateFilterOnePopulation}
        />
      )
    })
    nCasesFilterDisplay = {
      showFilter: true,
      filters: (
        <>
          <ColumnGroupIndividualFilters
            columns={columns.find(col => col.columnGroupName === ColumnGroupName.NCases).columns}
          />
          {nCasesPopulationFilterElems}
        </>
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
    const nControlsPopulationFilterElems = Object.entries(nControlsFilters).map(([population, filterValue]) => {
      return (
        <NumberRangeColumnFilter key={population}
          label={population}
          population={population as PopulationCode}
          globalMinValue={nControlsPerPopulationExtremums[population as PopulationCode].min}
          globalMaxValue={nControlsPerPopulationExtremums[population as PopulationCode].max}
          metric={RangeFilterMetric.NControls}
          filterValue={filterValue}
          disableFilter={disableFilterOnePopulation}
          updateFilter={updateFilterOnePopulation}
        />
      )
    })
    nControlsFilterDisplay = {
      showFilter: true,
      filters: (
        <>{nControlsPopulationFilterElems}</>
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
    const saigeHeritabilityPopulationFilterElems = Object.entries(saigeHeritabilityFilters).map(([population, filterValue]) => {
      return (
        <NumberRangeColumnFilter key={population}
          label={population}
          population={population as PopulationCode}
          globalMinValue={saigeHeritabilityPerPopulationExtremums[population as PopulationCode].min}
          globalMaxValue={saigeHeritabilityPerPopulationExtremums[population as PopulationCode].max}
          metric={RangeFilterMetric.SaigeHeritability}
          filterValue={filterValue}
          disableFilter={disableFilterOnePopulation}
          updateFilter={updateFilterOnePopulation}
        />
      )
    })
    saigeHeritabilityFilterDisplay = {
      showFilter: true,
      filters: (
        <>{saigeHeritabilityPopulationFilterElems}</>
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
    const lambdaGcPopulationFilterElems = Object.entries(lambdaGcFilters).map(([population, filterValue]) => {
      return (
        <NumberRangeColumnFilter key={population}
          label={population}
          population={population as PopulationCode}
          globalMinValue={lambdaGcPerPopulationExtremums[population as PopulationCode].min}
          globalMaxValue={lambdaGcPerPopulationExtremums[population as PopulationCode].max}
          metric={RangeFilterMetric.LambdaGc}
          filterValue={filterValue}
          disableFilter={disableFilterOnePopulation}
          updateFilter={updateFilterOnePopulation}
        />
      )
    })
    lambdaGcFilterDisplay = {
      showFilter: true,
      filters: (
        <> {lambdaGcPopulationFilterElems} </>
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
