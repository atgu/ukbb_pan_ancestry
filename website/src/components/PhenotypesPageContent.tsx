import React, { useCallback, useEffect, useMemo, useReducer, useRef, useState } from 'react'
import BrowserOnly from '@docusaurus/BrowserOnly';
import {  useTable, useFilters, useBlockLayout, useGlobalFilter  } from "react-table";
import {  SelectColumnFilter } from "../components/SelectColumnFilter";
import {  FixedSizeList } from "react-window";
import { Datum } from '../components/types';
import { scrollbarWidth } from '../components/scrollbarWidth';
import { Table, TableHead, TableRow, TableCell, TableBody, createMuiTheme, ThemeProvider, makeStyles, Theme, } from '@material-ui/core';
import {  DownloadLinkCell, downloadLinkCellMaxWidth } from "../components/DownloadLinkCell";
import {  CopyLinkCell, copyLinkCellMaxWidth } from "../components/CopyLinkCell";
import {  customIncludesFilterFn } from "../components/customIncludesFilterFn";
import { fuzzyTextFilterFunction, fuzzyTextGlobalFilterFunction } from '../components/fuzzyTextFilterFunction';
import { TextColumnFilter } from '../components/TextColumnFilter';

import data from "../data.json"
import { rangeFilterFunction } from '../components/rangeFilterFunction';
import phenotypesStyles from "./phenotypes.module.css"
import { TruncatedTextCell } from '../components/TruncatedTextCell';
import { PopulationCell } from '../components/PopulationCell';
import {  commonPopulations, PerPopulationMetrics, PopulationCode, } from "../components/populations";
import { NumberRangeColumnFilter } from '../components/NumberRangeColumnFilter';
import { PhenotypeFilters, PerPopulationExtremums } from '../components/PhenotypeFilters';
import { PopulationsFilter } from '../components/PopulationsFilter';
import { populationsFilterFunction } from '../components/populationsFilterFunction';
import {  BarChartCell, width as chartCellWidth, height as chartCellHeight } from "../components/BarChartCell";
import {  ScatterPlotCell} from "../components/ScatterPlotCell";
import {  ActionType, ColumnGroupName, initialState, PerPopulationMetricsVisibility, RangeFilterMetric, reducer } from "../components/phenotypesReducer";
import clsx from "clsx"
import { PopuplationHeader } from './PopulationHeader';
import { determineExtremums, maxSaigeHeritabilityValue } from './determineExtremums';
import { format } from 'd3-format';
import "regenerator-runtime/runtime";

const useStyles = makeStyles((theme: Theme) => ({
  oddTableRow: {
    // Need to double up on class names to override Infima's default styling:
    "&&": {
      backgroundColor: theme.palette.action.hover,
    }
  },
  tableCellNoPadding: {
    padding: "0",
  }
}))

const DefaultColumnFilter = () => null
const numCasesColumnWidth = 100
const decimalNumbersColumnWidth = 120

// Handle columns that can have "NA" values by replacing "NA" with `null`.
// Note: if `accessor` is a function, the `id` field must be present.
const getNaColumnProps = (fieldName) => ({
  id: fieldName,
  accessor: (row) => {
    const value = row[fieldName]
    return (value === "NA") ? null : value
  }
})

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


export const PhenotypesPageContent = () => {

  const [state, dispatch] = useReducer(reducer, initialState)
  const {
    columnGroupVisibilities, perPopulationMetricsVisibilities,
    n_cases: nCasesFilters,
    n_controls: nControlsFilters,
    saige_heritability: saigeHeritabilityFilters,
    lambda_gc: lambdaGcFilters,
  } = state;
  const classes = useStyles()

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
    data,
    perPopulationMetricsVisibilities,
    nCasesFilters,
    nControlsFilters,
    saigeHeritabilityFilters,
    lambdaGcFilters,
  }) , [data, perPopulationMetricsVisibilities, nCasesFilters, nControlsFilters, saigeHeritabilityFilters, lambdaGcFilters])

  const columns = useMemo(
    () => {
      const columnVisibilities = columnGroupVisibilities
      let columns = [
        {
          Header: "Description",
          accessor: "description",
          filter: fuzzyTextFilterFunction,
          Filter: TextColumnFilter,
          width: 300,
          Cell: TruncatedTextCell
        },
      ]
      if (columnVisibilities.analysis) {
        columns = [
          ...columns,
          {
            Header: "Analysis",
            columnGroupName: ColumnGroupName.Analysis,
            columns: [
              {
                Header: "Trait type", accessor: "trait_type",
                Filter: SelectColumnFilter,
                filter: customIncludesFilterFn,
                width: 120,
              },
              {
                Header: "Phenocode", accessor: "phenocode",
                filter: fuzzyTextFilterFunction,
                Filter: TextColumnFilter,
                Cell: TruncatedTextCell,
                width: 120,
              },
            ]
          },
        ]
      }
      if (columnVisibilities.description) {
        columns = [
          ...columns,
          {
            Header: "Description",
            columnGroupName: ColumnGroupName.Description,
            columns: [
              {
                Header: "Category",
                accessor: "category",
                filter: fuzzyTextFilterFunction,
                Filter: TextColumnFilter,
                width: 400,
                Cell: TruncatedTextCell
              },
              {
                Header: "Populations",
                accessor: "pops",
                filter: populationsFilterFunction,
                Filter: PopulationsFilter,
                Cell: PopulationCell,
                Header: PopuplationHeader,
              }
            ]
          },
        ]
      }
      if (columnVisibilities.nCases) {
        columns = [
          ...columns,
          {
            Header: "N Cases",
            columnGroupName: ColumnGroupName.NCases,
            columns: [
              {
                Header: "Both sexes",
                ...getNaColumnProps("n_cases_full_cohort_both_sexes"),
                Filter: NumberRangeColumnFilter,
                filter: rangeFilterFunction,
                width: numCasesColumnWidth,
              },
              {
                Header: "Females",
                ...getNaColumnProps("n_cases_full_cohort_females"),
                Filter: NumberRangeColumnFilter,
                filter: rangeFilterFunction,
                width: numCasesColumnWidth,
              },
              {
                Header: "Males",
                ...getNaColumnProps("n_cases_full_cohort_males"),
                Filter: NumberRangeColumnFilter,
                filter: rangeFilterFunction,
                width: numCasesColumnWidth,
              },
              {
                Header: "Per Population",
                ...getPerPopulationMetrics("n_cases", perPopulationMetricsVisibilities),
                Cell: BarChartCell,
                width: chartCellWidth,
                upperThreshold: globalNCasesUpperThreshold,
                labelFormatter: format(","),
              }
            ]
          }
        ]
      }
      if (columnVisibilities.nControls) {
        columns = [
          ...columns,
          {
            Header: "N Controls",
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
            Header: "Saige heritability",
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
            Header: "Lambda GC",
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
                Header: "tsv",
                accessor: "aws_link",
                Cell: DownloadLinkCell,
                width: downloadLinkCellMaxWidth,
                disableFilters: true,
              },
              {
                Header: "tbi",
                accessor: "aws_link_tabix",
                Cell: DownloadLinkCell,
                width: downloadLinkCellMaxWidth,
                disableFilters: true,
              },
              {
                Header: "wget tsv",
                accessor: "wget",
                Cell: CopyLinkCell,
                width: copyLinkCellMaxWidth,
                disableFilters: true,
              },
              {
                Header: "wget tbi",
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
              {Header: "tsv", accessor: "md5_hex", disableFilters: true, width: columnWidth},
              {Header: "tbi", accessor: "md5_hex_tabix", disableFilters: true, width: columnWidth},
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
        // {id: "n_cases_full_cohort_both_sexes", value: undefined},
        // {id: "n_cases_full_cohort_females", value: undefined},
        // {id: "n_cases_full_cohort_males", value: undefined},
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
    visibleColumns,
    state: reactTableState,
    columns: outputColumns,
  } = useTable<Datum>({
    columns,
    data: filteredByPopulationRangeFilters as any,
    initialState: initialReactTableState,
    defaultColumn,
    globalFilter: fuzzyTextGlobalFilterFunction,
  },
  useFilters,
  useGlobalFilter,
  useBlockLayout)

  const headerGroupElems = headerGroups.map(headerGroup => {
    const headerElems = headerGroup.headers.map(column => {
      return (
        <TableCell align="center" {...column.getHeaderProps()} className={classes.tableCellNoPadding}>
          {column.render("Header")}
        </TableCell>
      )
    })
    return (
      <TableRow {...headerGroup.getHeaderGroupProps()}>
        {headerElems}
      </TableRow>
    )
  })



  const RenderRow = useCallback(({index, style}) => {
    const row = rows[index]
    prepareRow(row)
    const cellElems = row.cells.map(cell => {
      return (
        <TableCell {...cell.getCellProps()} className={classes.tableCellNoPadding}>
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
  }, [prepareRow, rows])

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

  const scrollBarSizeRef = useRef<number | undefined>(undefined)
  useEffect(() => {
    if (scrollBarSizeRef.current === undefined) {
      scrollBarSizeRef.current = scrollbarWidth()
    }
  })
  return (
    <>
      <header>
        <div className={`container ${phenotypesStyles.titleContainer}`}>
          <h1 className="page-title">Phenotypes</h1>
        </div>
      </header>
      <main>
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
            />
          </div>
          <div style={{overflowX: "auto"}}>
            <Table size="small" {...getTableProps()}>
              <TableHead>
                {headerGroupElems}
              </TableHead>
              <BrowserOnly>
              {
                () => {
                  if (scrollBarSizeRef.current === undefined) {
                    return null
                  } else {
                    return (
                      <div {...getTableBodyProps()}>
                        <FixedSizeList
                          height={700}
                          itemCount={rows.length}
                          itemSize={chartCellHeight}
                          width={totalColumnsWidth + scrollBarSizeRef.current}
                          innerElementType={TableBody}
                        >
                          {RenderRow}
                        </FixedSizeList>
                      </div>
                    )
                  }
                }
              }
              </BrowserOnly>
            </Table>
          </div>
        </div>
      </main>
    </>
  )
}
