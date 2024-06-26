import React, { useCallback, useEffect, useMemo, useReducer, useState } from 'react'
import BrowserOnly from '@docusaurus/BrowserOnly';
import {  useTable, useFilters, useBlockLayout, useGlobalFilter} from "react-table";
import {  FixedSizeList } from "react-window";
import { Datum } from '../components/types';
import { scrollbarWidth } from '../components/scrollbarWidth';
import { Table, TableHead, TableRow, TableCell, TableBody, makeStyles, Theme, useTheme, useMediaQuery, } from '@material-ui/core';
import {  DownloadLinkCell, downloadLinkCellMaxWidth } from "../components/DownloadLinkCell";
import {  CopyLinkCell, copyLinkCellMaxWidth } from "../components/CopyLinkCell";
import { fuzzyTextFilterFunction, fuzzyTextGlobalFilterFunction } from '../components/fuzzyTextFilterFunction';
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
import { CenteredHeaderCell, SaigeHeritabilityHeaderCell } from './CenteredHeaderCell';
import { Option } from './DropdownFilter';
import Skeleton from '@material-ui/lab/Skeleton';
import { processPhenotypeDescription } from './descriptionAccessor';


const overallPageMaxWidth = "95%"
export const docusaurusLayoutWrapperClassName = "phenotypes-docusaurus-wrapper"
// How tall docusaurus's "Layout" element should be in mobile:
const mobileMinLayoutContainerHeight = 800; // in px

const useStyles = makeStyles((theme: Theme) => ({
  "@global": {
    [`.${docusaurusLayoutWrapperClassName}`]: {
      display: "flex",
      flexDirection: "column",
      [theme.breakpoints.down("sm")]: {
        flexBasis: `${mobileMinLayoutContainerHeight}px`,
      }
    }
  },
  oddTableRow: {
    // Need to double up on class names to override Infima's default styling:
    "&&": {
      backgroundColor: theme.palette.action.hover,
    }
  },
  main: {
    flexGrow: 1,
    display: "flex",
  },
  tableContainer: {
    overflowX: "auto",
    height: "100%",
    [theme.breakpoints.up("md")]: {
      height: "100%"
    },
    [theme.breakpoints.down("sm")]: {
      flexGrow: 1,
      marginTop: theme.spacing(1),
    }
  },
  container: {
    maxWidth: overallPageMaxWidth,
    "--ifm-table-stripe-background": "var(--ifm-table-background)",
    [theme.breakpoints.up("md")]: {
      display: "grid",
      gridTemplateColumns: "300px 1fr",
      columnGap: "20px",
    },
    [theme.breakpoints.down("sm")]: {
      display: "flex",
      flexDirection: "column",

    },
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

const getPhenotypeDescription = (row: Datum): Description => {
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

  const theme = useTheme()
  const isLargeScreen = useMediaQuery(theme.breakpoints.up("md"))
  const [state, dispatch] = useReducer(reducer, initialState)
  const [data, setData] = useState<Datum[] | undefined>(undefined)
  useEffect(() => {
    const fetchData = async () => {
      try {
        const loadedModule: {default: Datum[]} = await import( /* webpackChunkName: "data" */"../data.json")
        const withProcessedDescription = loadedModule.default.map(elem => processPhenotypeDescription(elem))
        setData(withProcessedDescription)
      } catch (e) {
        console.error(e)
      }
    }
    fetchData()
  }, [])
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

  type DerivedValues = undefined | (
    Omit<ReturnType<typeof determineExtremums>, "filteredByPopulationRangeFilters"> & { traitTypeFilterOptions: Option[] }
  )

  const {filteredData, derivedValues} = useMemo(() => {
    let derivedValues: DerivedValues
    let filteredData: Datum[]
    if (data === undefined) {
      filteredData = []
      derivedValues = undefined
    } else {
      const traitTypeValues = [...new Set(data.map(row => (row["trait_type"] as string)).filter(Boolean))] as string[]
      const traitTypeFilterOptions = traitTypeValues.map(val => ({label: val, value: val}))
      let filteredByDescriptionFilters = data
      if (sexFilterValue !== undefined) {
        filteredByDescriptionFilters = filteredByDescriptionFilters.filter(row => row["pheno_sex"] === sexFilterValue)
      }
      if (traitTypeFilterValue !== undefined) {
        filteredByDescriptionFilters = filteredByDescriptionFilters.filter(row => row["trait_type"] === traitTypeFilterValue)
      }
      if (descriptionFilterValue !== undefined) {
        filteredByDescriptionFilters = fuzzyTextFilterFunction(filteredByDescriptionFilters, "description", descriptionFilterValue)
      }
      const {filteredByPopulationRangeFilters, ...rest} = determineExtremums({
        data: filteredByDescriptionFilters,
        perPopulationMetricsVisibilities,
        nCasesFilters,
        nControlsFilters,
        saigeHeritabilityFilters,
        lambdaGcFilters,
      })
      derivedValues = {
        ...rest,
        traitTypeFilterOptions,
      }
      filteredData = filteredByPopulationRangeFilters
    }
    return {filteredData , derivedValues }
  }, [
    data, sexFilterValue, traitTypeFilterValue, descriptionFilterValue,
    perPopulationMetricsVisibilities, nCasesFilters, nControlsFilters, saigeHeritabilityFilters, lambdaGcFilters
  ])

  const columns = useMemo(
    () => {

      if (derivedValues === undefined) {
        return []
      } else {
        const {
          globalNCasesUpperThreshold, globalNControlsUpperThreshold,
          globalLambdaGc,
        } = derivedValues;
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
              Header: SaigeHeritabilityHeaderCell,
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
    }
  , [
    columnGroupVisibilities, perPopulationMetricsVisibilities,
    derivedValues,
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
    setGlobalFilter,
    state: reactTableState,
    columns: outputColumns,
    // These are the ones that are filtered by filters that are registered with react-table:
    rows: fullyFilteredRows,
  } = useTable<Datum>({
    columns,
    initialState: initialReactTableState,
    defaultColumn,
    globalFilter: fuzzyTextGlobalFilterFunction,
    // These are the ones that are filtered by filteres outside of react-table:
    data: filteredData,
  },
  useGlobalFilter,
  useFilters,
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
        <div className="container" style={{maxWidth: overallPageMaxWidth}}>
          <h1 className="page-title">Phenotypes</h1>
        </div>
      </header>
      <main className={classes.main}>
        <BrowserOnly>
        {
          () => {
            if (derivedValues === undefined) {
              return (
                <div className={clsx("container", classes.container)}>
                  <Skeleton variant="rect" height="70vh"/>
                  <Skeleton variant="rect" height="70vh"/>
                </div>
              )
            } else {
              const {
                perPopulationNCasesExtremums,
                perPopulationNControlsExtremums,
                perPopulationSaigeHeritabilityExtremums,
                perPopulationLambdaGcExtremums,
              } = derivedValues;
              return (
                    <div className={clsx("container", classes.container)}>
                      <div>
                        <PhenotypeFilters
                          isLargeScreen={isLargeScreen}
                          columnVisibilities={columnGroupVisibilities}
                          setColumnVisibilities={setColumnGroupVisibilites}
                          columns={outputColumns}
                          preGlobalFilteredRows={fullyFilteredRows}
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
                          traitTypeFilterOptions={derivedValues.traitTypeFilterOptions}
                          traitTypeFilterValue={traitTypeFilterValue}
                          setTraitTypeFilterValue={setTraitTypeFilterValue}
                          recordsCount={fullyFilteredRows.length}
                          descriptionFilterValue={descriptionFilterValue}
                          setDescriptionFilterValue={setDescriptionFilterValue}
                        />
                      </div>
                      <div className={classes.tableContainer}>
                        <AutoSizer>
                          {(dimensions) => {
                              const pageHeaderHeight = 106
                              const mobileTableFilterAndCollapsedControlHeight = 122
                              const mobileTableSeparationFromFooter = 20
                              let tableBody: React.ReactNode
                              const tableBodyHeight = isLargeScreen ?
                                dimensions.height - tableHeaderHeight - 50 :
                                // Note: subtract 10 to give some separation between table and footer:
                                mobileMinLayoutContainerHeight - pageHeaderHeight -
                                  mobileTableFilterAndCollapsedControlHeight - tableHeaderHeight - mobileTableSeparationFromFooter
                              if (scrollBarSize === undefined) {
                                tableBody = null
                              } else {
                                tableBody = (
                                  <div {...getTableBodyProps()}>
                                    <FixedSizeList
                                      // Subtract away some extra space for safety to minimize the chance of vertical scroll bar being added:
                                      height={tableBodyHeight}
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
        }
      </BrowserOnly>
      </main>
    </>
  )
}
