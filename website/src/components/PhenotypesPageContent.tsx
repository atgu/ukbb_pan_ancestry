import React, { useCallback, useEffect, useMemo, useReducer, useState } from 'react'
import BrowserOnly from '@docusaurus/BrowserOnly';
import {  useTable, useFilters, useBlockLayout, useGlobalFilter} from "react-table";
import {  FixedSizeList } from "react-window";
import { Datum } from '../components/types';
import { scrollbarWidth } from '../components/scrollbarWidth';
import { Table, TableHead, TableRow, TableCell, TableBody, makeStyles, Theme, } from '@material-ui/core';
import {  DownloadLinkCell, downloadLinkCellMaxWidth } from "../components/DownloadLinkCell";
import {  CopyLinkCell, copyLinkCellMaxWidth } from "../components/CopyLinkCell";
import { fuzzyTextFilterFunction, fuzzyTextGlobalFilterFunction } from '../components/fuzzyTextFilterFunction';
import data from "../data.json"
import phenotypesStyles from "./phenotypes.module.css"
import { PopulationCell } from '../components/PopulationCell';
import {  commonPopulations, PerPopulationMetrics, PopulationCode, } from "../components/populations";
import { PhenotypeFilters} from '../components/PhenotypeFilters';
import { PopulationsFilter } from '../components/PopulationsFilter';
import { populationsFilterFunction } from '../components/populationsFilterFunction';
import {  BarChartCell, width as chartCellWidth, height as chartCellHeight } from "../components/BarChartCell";
import {  ScatterPlotCell} from "../components/ScatterPlotCell";
import {  ActionType, ColumnGroupName, initialState, PerPopulationMetricsVisibility, RangeFilterMetric, reducer } from "../components/phenotypesReducer";
import clsx from "clsx"
import { PopuplationHeader } from './PopulationHeader';
import { determineExtremums, maxSaigeHeritabilityValue } from './determineExtremums';
import { format } from 'd3-format';
import { Description, DescriptionCell, width as descriptionCellWidth } from './DescriptionCell';
import {  AutoSizer } from "react-virtualized";
import { CenteredHeaderCell } from './CenteredHeaderCell';


const useStyles = makeStyles((theme: Theme) => ({
  oddTableRow: {
    // Need to double up on class names to override Infima's default styling:
    "&&": {
      backgroundColor: theme.palette.action.hover,
    }
  },
  main: {
    flexGrow: 1,
    display: "flex"
  },
  tableContainer: {
    overflowX: "auto",
    height: "100%"
  }
}))

const DefaultColumnFilter = () => null

const getPerPopulationMetrics = (
    metric_prefix: string,
    populationMetricsVisibilities: PerPopulationMetricsVisibility,
  ) => {
  return {
    id: `${metric_prefix}_per_population`,
    accessor: (row) => {
      const result: PerPopulationMetrics = new Map()
      for (const populationCode of commonPopulations) {
        const key = `${metric_prefix}_${populationCode}`
        const populationValue = row[key]
        if (typeof populationValue === "number" && populationMetricsVisibilities[populationCode] === true) {
          result.set(populationCode, populationValue)
        }
      }
      return result
    }
  }
}

const getPhenotypeDescription = (row): Description => {
  const sex = (row.pheno_sex) ? row.pheno_sex : undefined
  const code = (row.phenocode) ? row.phenocode : undefined
  const traitType = (row.trait_type) ? row.trait_type : undefined
  return {
    name: row.description,
    category: row.category,
    sex,
    code,
    traitType,
  }
}


export const PhenotypesPageContent = () => {

  const [state, dispatch] = useReducer(reducer, initialState)
  const {
    columnGroupVisibilities, perPopulationMetricsVisibilities,
    n_cases: nCasesFilters,
    n_controls: nControlsFilters,
    saige_heritability: saigeHeritabilityFilters,
    lambda_gc: lambdaGcFilters,
    sexFilterValue, traitTypeFilterValue,
    descriptionFilterValue,
  } = state;
  const classes = useStyles()

  const {filteredByDescriptionFilters, traitTypeFilterOptions} = useMemo(() => {
    const traitTypeValues = [...new Set(data.map(row => (row["trait_type"] as string)).filter(Boolean))] as string[]
    const traitTypeFilterOptions = traitTypeValues.map(val => ({label: val, value: val}))
    let filtered = data
    if (sexFilterValue !== undefined) {
      filtered = filtered.filter(row => row["pheno_sex"] === sexFilterValue)
    }
    if (traitTypeFilterValue !== undefined) {
      filtered = filtered.filter(row => row["trait_type"] === traitTypeFilterValue)
    }
    if (descriptionFilterValue !== undefined) {
      filtered = fuzzyTextFilterFunction(filtered, "description", descriptionFilterValue)
    }
    return {
      filteredByDescriptionFilters: filtered,
      traitTypeFilterOptions,
    }
  }, [sexFilterValue, traitTypeFilterValue, descriptionFilterValue])

  const {
    globalNCasesUpperThreshold,
    perPopulationNCasesExtremums,
    globalNControlsUpperThreshold,
    perPopulationNControlsExtremums,
    perPopulationSaigeHeritabilityExtremums,
    globalLambdaGc,
    filteredByPopulationRangeFilters,
    perPopulationLambdaGcExtremums,
  } = useMemo( () => determineExtremums({
    data: filteredByDescriptionFilters,
    perPopulationMetricsVisibilities,
    nCasesFilters,
    nControlsFilters,
    saigeHeritabilityFilters,
    lambdaGcFilters,
  }) , [filteredByDescriptionFilters, perPopulationMetricsVisibilities, nCasesFilters, nControlsFilters, saigeHeritabilityFilters, lambdaGcFilters])

  const columns = useMemo(
    () => {
      const columnVisibilities = columnGroupVisibilities
      let columns = [
        {
          Header: CenteredHeaderCell,
          width: descriptionCellWidth,
          accessor: getPhenotypeDescription,
          id: "description",
          columnTitle: "Description",
          Cell: DescriptionCell
        },
      ]
      if (columnVisibilities.populations) {
        columns = [
          ...columns,
          {
            columnGroupName: ColumnGroupName.Populations,
            accessor: "pops",
            filter: populationsFilterFunction,
            Filter: PopulationsFilter,
            Cell: PopulationCell,
            Header: PopuplationHeader,
          },
        ]
      }
      if (columnVisibilities.nCases) {
        columns = [
          ...columns,
          {
            Header: CenteredHeaderCell,
            columnTitle: "N Cases",
            ...getPerPopulationMetrics("n_cases", perPopulationMetricsVisibilities),
            Cell: BarChartCell,
            width: chartCellWidth,
            upperThreshold: globalNCasesUpperThreshold,
            labelFormatter: format(","),
            columnGroupName: ColumnGroupName.NCases,
          }
        ]
      }
      if (columnVisibilities.nControls) {
        columns = [
          ...columns,
          {
            Header: CenteredHeaderCell,
            columnTitle: "N Controls",
            columnGroupName: ColumnGroupName.NControls,
            ...getPerPopulationMetrics("n_controls", perPopulationMetricsVisibilities),
            Cell: BarChartCell,
            width: chartCellWidth,
            upperThreshold: globalNControlsUpperThreshold,
            labelFormatter: format(","),
          },
        ]
      }
      if (columnVisibilities.saigeHeritability) {
        columns = [
          ...columns,
          {
            Header: CenteredHeaderCell,
            columnTitle: "Saige heritability",
            columnGroupName: ColumnGroupName.SaigeHeritability,
            ...getPerPopulationMetrics("saige_heritability", perPopulationMetricsVisibilities),
            Cell: BarChartCell,
            width: chartCellWidth,
            upperThreshold: maxSaigeHeritabilityValue,
            labelFormatter: format(".3f"),
          }
        ]
      }
      if (columnVisibilities.lambdaGc) {
        columns = [
          ...columns,
          {
            Header: CenteredHeaderCell,
            columnTitle: "Lambda GC",
            columnGroupName: ColumnGroupName.LambdaGc,
            ...getPerPopulationMetrics("lambda_gc", perPopulationMetricsVisibilities),
            Cell: ScatterPlotCell,
            width: chartCellWidth,
            minValue: globalLambdaGc.min,
            maxValue: globalLambdaGc.max,
            xAxisConfig: {isShown: true, yIntercept: 1},
            labelFormatter: format(".3f"),
          }
        ]
      }
      if (columnVisibilities.downloads) {
        columns = [
          ...columns,
          {
            Header: "Downloads",
            columnGroupName: ColumnGroupName.Downloads,
            columns: [
              {
                Header: CenteredHeaderCell,
                columnTitle: "tsv",
                accessor: "aws_link",
                Cell: DownloadLinkCell,
                width: downloadLinkCellMaxWidth,
                disableFilters: true,
              },
              {
                Header: CenteredHeaderCell,
                columnTitle: "tbi",
                accessor: "aws_link_tabix",
                Cell: DownloadLinkCell,
                width: downloadLinkCellMaxWidth,
                disableFilters: true,
              },
              {
                Header: CenteredHeaderCell,
                columnTitle: "wget tsv",
                accessor: "wget",
                Cell: CopyLinkCell,
                width: copyLinkCellMaxWidth,
                disableFilters: true,
              },
              {
                Header: CenteredHeaderCell,
                columnTitle: "wget tbi",
                accessor: "wget_tabix",
                Cell: CopyLinkCell,
                width: copyLinkCellMaxWidth,
                disableFilters: true,
              },
            ]
          }
        ]
      }
      if (columnVisibilities.md5) {
        const columnWidth = 320;
        columns = [
          ...columns,
          {
            Header: "MD5",
            columnGroupName: ColumnGroupName.Md5,
            columns: [
              {
                Header: CenteredHeaderCell,
                columnTitle: "tsv",
                accessor: "md5_hex",
                disableFilters: true,
                width: columnWidth
              },
              {
                Header: CenteredHeaderCell,
                columnTitle: "tbi",
                accessor: "md5_hex_tabix",
                disableFilters: true,
                width: columnWidth
              },
            ]
          }
        ]
      }
      return columns
    }
  , [
    columnGroupVisibilities, perPopulationMetricsVisibilities,
    globalNCasesUpperThreshold, globalNControlsUpperThreshold,
    globalLambdaGc.min, globalLambdaGc.max,
  ])
  const initialReactTableState = useMemo(() => {
    return {
      globalFilter: "",
      filters: [
        {id: "trait_type", value: ""},
        {id: "pheno_sex", value: ""},
      ]
    }
  }, [])

  const defaultColumn = React.useMemo(
    () => ({
      Filter: DefaultColumnFilter,
    }),
    []
  )
  const {
    getTableProps, getTableBodyProps, headerGroups, rows, prepareRow,
    totalColumnsWidth,
    preGlobalFilteredRows,
    setGlobalFilter,
    state: reactTableState,
    columns: outputColumns,
  } = useTable<Datum>({
    columns,
    data: filteredByPopulationRangeFilters,
    initialState: initialReactTableState,
    defaultColumn,
    globalFilter: fuzzyTextGlobalFilterFunction,
  },
  useFilters,
  useGlobalFilter,
  useBlockLayout)

  const headerGroupElems = headerGroups.map((headerGroup, headerGroupIndex) => {
    const headerElems = headerGroup.headers.map((column, columnIndex) => {
      return (
        <TableCell align="center" {...column.getHeaderProps()} key={columnIndex}>
          {column.render("Header")}
        </TableCell>
      )
    })
    return (
      <TableRow {...headerGroup.getHeaderGroupProps()} key={headerGroupIndex}>
        {headerElems}
      </TableRow>
    )
  })



  const RenderRow = useCallback(({index, style}) => {
    const row = rows[index]
    prepareRow(row)
    const cellElems = row.cells.map((cell, cellIndex) => {
      return (
        <TableCell {...cell.getCellProps()} key={cellIndex}>
          {cell.render("Cell")}
        </TableCell>
      )
    })
    return (
      <TableRow
        className={clsx({
          [classes.oddTableRow]: index % 2 == 1
        })}
        {...row.getRowProps({style})}
      >
        {cellElems}
      </TableRow>
    )
  }, [prepareRow, rows, classes.oddTableRow])

  const setColumnGroupVisibilites = useCallback((columnGroupName: ColumnGroupName, isVisible: boolean) => dispatch({
    type: ActionType.SET_COLUMN_GROUP_VISIBILITY,
    payload: {
      columnGroup: columnGroupName,
      isVisible,
    }
  }), [])

  const setPerPoulationMetricsVisibilities = useCallback((population: PopulationCode, isVisible: boolean) => dispatch({
    type: ActionType.SET_POPULATION_METRICS_VISIBILITY,
    payload: { population, isVisible }
  }), [])

  const disableRangeFilterOnePopulation = useCallback((args: {metric: RangeFilterMetric, population: PopulationCode}) => dispatch({
    type: ActionType.DISABLE_FILTER_ONE_POPULATION,
    payload: {metric: args.metric, population: args.population}
  }), [])
  const updateRangeFilterOnePopulation = useCallback((args: {metric: RangeFilterMetric, population: PopulationCode, min: number, max: number}) => dispatch({
    type: ActionType.UPDATE_FILTER_ONE_POPULATION,
    payload: {
      metric: args.metric, population: args.population,
      min: args.min, max: args.max
    }
  }), [])

  const [scrollBarSize, setScrollBarSize] = useState<number | undefined>(undefined)
  useEffect(() => {
    if (scrollBarSize === undefined) {
      setScrollBarSize(scrollbarWidth)
    }
  }, [scrollBarSize])

  const setSexFilterValue = useCallback((value: string | undefined) => dispatch({
    type: ActionType.UPDATE_SEX_FILTER,
    payload: {value}
  }), [])
  const setTraitTypeFilterValue = useCallback((value: string | undefined) => dispatch({
    type: ActionType.UPDATE_TRAIT_TYPE_FILTER,
    payload: {value}
  }), [])
  const setDescriptionFilterValue = useCallback((value: string | undefined) => dispatch({
    type: ActionType.UPDATE_DESCRIPTION_FILTER,
    payload: {value}
  }), [])

  // From empirical measurement:
  const tableHeaderHeight = 54

  return (
    <>
      <header>
        <div className={`container ${phenotypesStyles.titleContainer}`}>
          <h1 className="page-title">Phenotypes</h1>
        </div>
      </header>
      <main className={classes.main}>
        <BrowserOnly>
        {
          () => {
            return (
                  <div className={`container ${phenotypesStyles.container}`}>
                    <div>
                      <PhenotypeFilters
                        columnVisibilities={columnGroupVisibilities}
                        setColumnVisibilities={setColumnGroupVisibilites}
                        columns={outputColumns}
                        preGlobalFilteredRows={preGlobalFilteredRows}
                        setGlobalFilter={setGlobalFilter}
                        globalFilter={reactTableState.globalFilter}
                        populationMetricsVisibilities={perPopulationMetricsVisibilities}
                        setPopulationMetricsVisibilities={setPerPoulationMetricsVisibilities}
                        nCasesFilters={nCasesFilters}
                        nCasesPerPopulationExtremums={perPopulationNCasesExtremums}
                        nControlsFilters={nControlsFilters}
                        nControlsPerPopulationExtremums={perPopulationNControlsExtremums}
                        saigeHeritabilityFilters={saigeHeritabilityFilters}
                        saigeHeritabilityPerPopulationExtremums={perPopulationSaigeHeritabilityExtremums}
                        lambdaGcFilters={lambdaGcFilters}
                        lambdaGcPerPopulationExtremums={perPopulationLambdaGcExtremums}
                        disableFilterOnePopulation={disableRangeFilterOnePopulation}
                        updateFilterOnePopulation={updateRangeFilterOnePopulation}
                        sexFilterValue={sexFilterValue}
                        setSexFilterValue={setSexFilterValue}
                        traitTypeFilterOptions={traitTypeFilterOptions}
                        traitTypeFilterValue={traitTypeFilterValue}
                        setTraitTypeFilterValue={setTraitTypeFilterValue}
                        recordsCount={filteredByPopulationRangeFilters.length}
                        descriptionFilterValue={descriptionFilterValue}
                        setDescriptionFilterValue={setDescriptionFilterValue}
                      />
                    </div>
                    <div className={classes.tableContainer}>
                      <AutoSizer>
                        {(dimensions) => {
                            let tableBody: React.ReactNode
                            if (scrollBarSize === undefined) {
                              tableBody = null
                            } else {
                              tableBody = (
                                <div {...getTableBodyProps()}>
                                  <FixedSizeList
                                    // Subtract away some extra space for safety to minimize the chance of vertical scroll bar being added:
                                    height={dimensions.height - tableHeaderHeight - 50}
                                    itemCount={rows.length}
                                    itemSize={chartCellHeight}
                                    width={totalColumnsWidth + scrollBarSize}
                                    innerElementType={TableBody}
                                  >
                                    {RenderRow}
                                  </FixedSizeList>
                                </div>
                              )
                            }
                          return (
                            <Table {...getTableProps()} >
                              <TableHead>
                                {headerGroupElems}
                              </TableHead>
                              {tableBody}

                            </Table>
                          )
                        }}
                      </AutoSizer>
                    </div>
                  </div>
            )
          }
        }
      </BrowserOnly>
      </main>
    </>
  )
}
