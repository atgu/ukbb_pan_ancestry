import React, { useState } from "react"
import { Accordion, AccordionDetails, AccordionSummary, Card, CardContent, FormControlLabel, FormGroup, makeStyles, Paper, Switch, Theme, Typography} from "@material-ui/core"
import { useCallback } from "react"
import { ColumnInstance, UseGlobalFiltersInstanceProps } from "react-table"
import { Datum } from "./types"
import { ExpandMore } from "@material-ui/icons"
import { ColumnGroupIndividualFilters } from "./ColumnGroupIndividualFilters"
import {  ColumnGroupFilterGroup, FilterDisplay, preventEventPropagation } from "./ColumnGroupFilterGroup";
import { GlobalFilter } from "./GlobalFilter"
import { commonPopulations, PopulationCode, populationColorMapping } from "./populations"
import { PopulationSwitch } from "./PopulationSwitch"
import { ColumnGroupVisibility, PerPopulationMetricsVisibility, ColumnGroupName, PerPopulationRangeFilter, RangeFilterMetric} from "./phenotypesReducer"
import {  NumberRangeColumnFilter } from "./NumberRangeColumnFilter";
import { DropdownFilter, Option } from "./DropdownFilter"
import { BaseTextFilter } from "./BaseTextFilter"
import {  PopulationColorLegend } from "./PopulationColorLegend";


const useStyles = makeStyles((theme: Theme) => ({
  accordionTitleNoSwitch: {
    marginLeft: theme.spacing(6)
  },
  populationMetricsFilterGroup: {
    display: "grid",
    gridTemplateColumns: "repeat(3, 1fr)",
  },
  descriptionFormGroup: {
    width: "100%",
  },
  traitTypeFilter: {
    marginTop: theme.spacing(1),
  },
  columnGroupCard: {
    marginBottom: theme.spacing(1)
  }
}))

export type PerPopulationExtremums = {
  [K in PopulationCode]: {min: number, max: number}
}

enum AccordionName {
  Populations,
  Description,
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
  globalFilter: string | undefined
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
  setSexFilterValue: (value: string | undefined) => void
  sexFilterValue: string | undefined
  traitTypeFilterOptions: Option[]
  traitTypeFilterValue: string | undefined
  setTraitTypeFilterValue: (value: string | undefined) => void
  descriptionFilterValue: string | undefined
  setDescriptionFilterValue: (value: string | undefined) => void
  recordsCount: number
  isLargeScreen: boolean
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
    sexFilterValue, setSexFilterValue,
    traitTypeFilterOptions, traitTypeFilterValue, setTraitTypeFilterValue,
    descriptionFilterValue, setDescriptionFilterValue,
    recordsCount,
    isLargeScreen
  } = props;
  const classes = useStyles()

  const [expandedAccordion, setExpandedAccordion] = useState<AccordionName | undefined>(undefined)

  const handleColumnVisibilityChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    setColumnVisibilities(event.target.name as ColumnGroupName, event.target.checked)
  }, [setColumnVisibilities])

  const getAccordionChangeHandler = (accordionName: AccordionName) => (_: unknown, isExpanded: boolean) => setExpandedAccordion(isExpanded ? accordionName : undefined)

  const descriptionFilter = (
    <FormGroup className={classes.descriptionFormGroup}>
      <DropdownFilter
        options={[
          {value: "both_sexes", label: "both"},
          {value: "females", label: "female"},
          {value: "males", label: "male"},
        ]}
        setFilterValue={setSexFilterValue}
        filterValue={sexFilterValue}
        label="Sex"
      />
      <DropdownFilter
        className={classes.traitTypeFilter}
        options={traitTypeFilterOptions}
        setFilterValue={setTraitTypeFilterValue}
        filterValue={traitTypeFilterValue}
        label="Trait Type"
      />
      <BaseTextFilter
        filterValue={descriptionFilterValue}
        setFilterValue={setDescriptionFilterValue}
        label={`Search ${recordsCount} records`}
      />
    </FormGroup>

  )

  const descriptionAccordion = (
    <Accordion expanded={expandedAccordion === AccordionName.Description} onChange={getAccordionChangeHandler(AccordionName.Description)}>
      <AccordionSummary expandIcon={<ExpandMore/>} className={classes.accordionTitleNoSwitch}>
        Description
      </AccordionSummary>
      <AccordionDetails>
        {descriptionFilter}
      </AccordionDetails>
    </Accordion>
  )

  const downloadsFilter = (
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
  )

  const downloadsAccordion = (
    <Paper>
      <AccordionSummary>
        {downloadsFilter}
      </AccordionSummary>
    </Paper>
  )

  let populationsAccordionDetails: React.ReactNode
  if (columnVisibilities[ColumnGroupName.Populations]) {
    populationsAccordionDetails =
      columns.find(col => col.columnGroupName === ColumnGroupName.Populations)!.render("Filter")
  } else {
    populationsAccordionDetails = null
  }

  const populationsFilter = (
    <Accordion
      expanded={expandedAccordion === AccordionName.Populations}
      onChange={getAccordionChangeHandler(AccordionName.Populations)}>
      <AccordionSummary expandIcon={<ExpandMore/>}>
        <FormControlLabel
          onClick={preventEventPropagation}
          onFocus={preventEventPropagation}
          control={
            <Switch
              checked={columnVisibilities[ColumnGroupName.Populations]} name={ColumnGroupName.Populations}
              onChange={handleColumnVisibilityChange}
            />
          }
          label="Populations"
        />
      </AccordionSummary>
      <AccordionDetails>
        {populationsAccordionDetails}
      </AccordionDetails>
    </Accordion>
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

  const perPopulationMetricsVisibilitiesFilter = (
    <FormGroup className={classes.populationMetricsFilterGroup}>
      {populationMetricsVisibilitiesFilterElems}
    </FormGroup>
  )

  const populationMetricsVibilitiesAccordion = (
    <Accordion expanded={expandedAccordion === AccordionName.PerPopulationMetrics} onChange={getAccordionChangeHandler(AccordionName.PerPopulationMetrics)}>
      <AccordionSummary expandIcon={<ExpandMore/>} className={classes.accordionTitleNoSwitch}>
        Per-Population Metrics
      </AccordionSummary>
      <AccordionDetails>
        {perPopulationMetricsVisibilitiesFilter}
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
      filterElems: (
        <>
          <ColumnGroupIndividualFilters
            columns={[columns.find(col => col.columnGroupName === ColumnGroupName.NCases)]}
          />
          {nCasesPopulationFilterElems}
        </>
      )
    }
  } else {
    nCasesFilterDisplay = {showFilter: false}
  }
  const nCasesVisibilityControl = (
    <Switch
      checked={columnVisibilities[ColumnGroupName.NCases]} name={ColumnGroupName.NCases}
      onChange={handleColumnVisibilityChange}
    />

  )
  const nCasesAccordion = (
    <ColumnGroupFilterGroup
      visibilityControl={nCasesVisibilityControl }
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
      filterElems: (
        <>{nControlsPopulationFilterElems}</>
      )
    }
  } else {
    nControlsFilterDisplay = {showFilter: false}
  }
  const nControlsVisibilityControl = (
    <Switch
      checked={columnVisibilities[ColumnGroupName.NControls]} name={ColumnGroupName.NControls}
      onChange={handleColumnVisibilityChange}
    />
  )
  const nControlsAccordion = (
    <ColumnGroupFilterGroup
      visibilityControl={nControlsVisibilityControl }
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
      filterElems: (
        <>{saigeHeritabilityPopulationFilterElems}</>
      )
    }

  } else {
    saigeHeritabilityFilterDisplay = {showFilter: false}
  }
  const saigeHeritabilityVisibilityControl = (
    <Switch
      checked={columnVisibilities[ColumnGroupName.SaigeHeritability]} name={ColumnGroupName.SaigeHeritability}
      onChange={handleColumnVisibilityChange}
    />
  )
  const saigeHeritabilityAccordion = (
    <ColumnGroupFilterGroup
      visibilityControl={saigeHeritabilityVisibilityControl }
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
      filterElems: (
        <> {lambdaGcPopulationFilterElems} </>
      )
    }

  } else {
    lambdaGcFilterDisplay = {showFilter: false}
  }
  const lambdaGcVisibilityControl = (
    <Switch
      checked={columnVisibilities[ColumnGroupName.LambdaGc]} name={ColumnGroupName.LambdaGc}
      onChange={handleColumnVisibilityChange}
    />

  )
  const lambdaGcAccordion = (
    <ColumnGroupFilterGroup
      visibilityControl={lambdaGcVisibilityControl }
      label="Lambda GC"
      filterDisplay={lambdaGcFilterDisplay}
      isAccordionExpanded={expandedAccordion === AccordionName.LambdaGc}
      onAccordionChange={getAccordionChangeHandler(AccordionName.LambdaGc)}
    />
  )
  const md5VisibilityControl = (
    <Switch
      checked={columnVisibilities[ColumnGroupName.Md5]} name={ColumnGroupName.Md5}
      onChange={handleColumnVisibilityChange}
    />
  )
  const md5FilterPanel = (
    <Paper>
      <AccordionSummary>
        <FormControlLabel
          onClick={preventEventPropagation}
          onFocus={preventEventPropagation}
          control={md5VisibilityControl }
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

  if (isLargeScreen) {
    return (
      <div>
        {globalFilterElem}
        <Paper>
          <AccordionSummary className={classes.accordionTitleNoSwitch}>
            Population Colors
          </AccordionSummary>
          <AccordionDetails>
            <PopulationColorLegend labelFontSize={undefined} iconFontSize={undefined} gridJustifyItems={undefined}/>
          </AccordionDetails>
        </Paper>
        {descriptionAccordion}
        {populationsFilter}
        {downloadsAccordion}
        {populationMetricsVibilitiesAccordion}
        {nCasesAccordion}
        {nControlsAccordion}
        {saigeHeritabilityAccordion}
        {lambdaGcAccordion}
        {md5FilterPanel}
      </div>
    )
  } else {
    return (
      <div>
        {globalFilterElem}
        <Accordion >
          <AccordionSummary expandIcon={<ExpandMore/>}>
            Table Controls
          </AccordionSummary>
          <AccordionDetails style={{flexDirection: "column"}}>
            <Card className={classes.columnGroupCard}>
              <CardContent>
              <Typography >Population Colors</Typography>
              <PopulationColorLegend labelFontSize={undefined} iconFontSize={undefined} gridJustifyItems="flex-start"/>
              </CardContent>
            </Card>

            <Card className={classes.columnGroupCard}>
              <CardContent>
                <Typography>Description</Typography>
                {descriptionFilter}
              </CardContent>
            </Card>

            <Card className={classes.columnGroupCard}>
              <CardContent>
              <FormControlLabel
                control={
                  <Switch
                    checked={columnVisibilities[ColumnGroupName.Populations]} name={ColumnGroupName.Populations}
                    onChange={handleColumnVisibilityChange}
                  />
                }
                label="Population"
              />
              {populationsAccordionDetails}
              </CardContent>
            </Card>

            <Card className={classes.columnGroupCard}>
              <CardContent>
                {downloadsFilter}
              </CardContent>
            </Card>

            <Card className={classes.columnGroupCard}>
              <CardContent>
                <Typography>Per-Population Metrics Visibility</Typography>
                {perPopulationMetricsVisibilitiesFilter}
              </CardContent>
            </Card>

            <Card className={classes.columnGroupCard}>
              <CardContent>
                <FormControlLabel control={ nCasesVisibilityControl } label="N Cases Column"/>
                {nCasesFilterDisplay.showFilter ? nCasesFilterDisplay.filterElems : null}
              </CardContent>
            </Card>

            <Card className={classes.columnGroupCard}>
              <CardContent>
                <FormControlLabel control={ nControlsVisibilityControl } label={
                  <Typography>N Controls Column</Typography>
                }/>
                {nControlsFilterDisplay.showFilter ? nControlsFilterDisplay.filterElems : null}
              </CardContent>
            </Card>

            <Card className={classes.columnGroupCard}>
              <CardContent>
                <FormControlLabel control={ saigeHeritabilityVisibilityControl} label="Saige Heritability Column" />
                {saigeHeritabilityFilterDisplay.showFilter ? saigeHeritabilityFilterDisplay.filterElems : null}
              </CardContent>
            </Card>

            <Card className={classes.columnGroupCard}>
              <CardContent>
                <FormControlLabel control={ lambdaGcVisibilityControl} label="Lambda GC Column" />
                {lambdaGcFilterDisplay.showFilter ? lambdaGcFilterDisplay.filterElems : null}
              </CardContent>
            </Card>

            <Card className={classes.columnGroupCard}>
              <CardContent>
                <FormControlLabel control={ md5VisibilityControl} label="MD5 Column" />
              </CardContent>
            </Card>

          </AccordionDetails>
        </Accordion>
      </div>
    )
  }

}
