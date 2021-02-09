import React, { useCallback, useMemo, useState } from 'react'
import Layout from '@theme/Layout'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import {  useTable, useFilters, useBlockLayout, useGlobalFilter  } from "react-table";
import {  SelectColumnFilter } from "../components/SelectColumnFilter";
import {  FixedSizeList } from "react-window";
import { Datum } from '../components/types';
import { scrollbarWidth } from '../components/scrollbarWidth';
import { Table, TableHead, TableRow, TableCell, TableBody, } from '@material-ui/core';
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
import { GlobalFilter } from '../components/GlobalFilter';
import { NumberRangeColumnFilter } from '../components/NumberRangeColumnFilter';
import { ColumnGroupVisibility, PhenotypeFilters, PopulationVisibility } from '../components/PhenotypeFilters';
import { PopulationsFilter } from '../components/PopulationsFilter';
import { populationsFilterFunction } from '../components/populationsFilterFunction';
import {  SaigeHeritabilityCell, width as saigeHeritabilityWidth } from "../components/SaigeHeritabilityCell";
import min from "lodash/min"
import max from "lodash/max"

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
    populationMetricsVisibilities: PopulationVisibility,
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

const Phenotypes = () => {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context

  const [columnVisibilities, setColumnVisibilities] = useState<ColumnGroupVisibility>({
    downloads: false,
    description: false,
    nCases: true,
    nControls: false,
    saigeHeritability: false,
    lambdaGc: false,
    md5: false,
  })

  const [populationMetricsVisibilities, setPopulationMetricsVisibilities] = useState<PopulationVisibility>({
    [PopulationCode.AFR]: true,
    [PopulationCode.AMR]: true,
    [PopulationCode.CSA]: true,
    [PopulationCode.EAS]: true,
    [PopulationCode.EUR]: true,
    [PopulationCode.MID]: true,
  })

  const minSaigeHeritabilityValue = 0
  const maxSaigeHeritabilityValue = 1

  const {
    minPopulationNCasesValue, maxPopulationNCasesValue
  } = useMemo(() => {
    const allPopulationNCasesValues = []
    for (const datum of data) {
      for (const population of commonPopulations) {
        const populationNCasesValue = datum[`n_cases_${population}`]
        if (typeof populationNCasesValue === "number") {
          allPopulationNCasesValues.push(populationNCasesValue)
        }
      }
    }
    return {
      minPopulationNCasesValue: min(allPopulationNCasesValues),
      maxPopulationNCasesValue: max(allPopulationNCasesValues)
    }

  }, [data])

  const columns = useMemo(
    () => {
      let columns = [
        {
          Header: "Description",
          accessor: "description",
          filter: fuzzyTextFilterFunction,
          Filter: TextColumnFilter,
          width: 300,
          Cell: TruncatedTextCell
        },
        {
          Header: "Analysis",
          columnGroupVisibilityAttribName: "analysis",
          columns: [
            {
              Header: "Coding", accessor: "coding",
              width: 100,
            },
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
            {
              Header: "Sex", accessor: "pheno_sex",
              Filter: SelectColumnFilter,
              filter: customIncludesFilterFn,
              width: 120,
            },
            {
              Header: "Modifier", accessor: "modifier",
              disableFilters: true,
              width: 100,
            },
          ]
        },
      ]
      if (columnVisibilities.description) {
        columns = [
          ...columns,
          {
            Header: "Description",
            columnGroupVisibilityAttribName: "description",
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
            columnGroupVisibilityAttribName: "nCases",
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
                ...getPerPopulationMetrics("n_cases", populationMetricsVisibilities),
                Cell: SaigeHeritabilityCell,
                width: saigeHeritabilityWidth,
                materialUiNoPadding: true,
                minValue: minPopulationNCasesValue,
                maxValue: maxPopulationNCasesValue,
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
            columnGroupVisibilityAttribName: "nControls",
            columns: commonPopulations.filter(pop => populationMetricsVisibilities[pop] === true).map(pop => ({
              Header: pop,
              ...getNaColumnProps(`n_controls_${pop}`),
              Filter: NumberRangeColumnFilter,
              filter: rangeFilterFunction,
              width: numCasesColumnWidth,
            }))
          }
        ]
      }
      if (columnVisibilities.saigeHeritability) {
        columns = [
          ...columns,
          {
            Header: "Saige heritability",
            columnGroupVisibilityAttribName: "saigeHeritability",
            ...getPerPopulationMetrics("saige_heritability", populationMetricsVisibilities),
            Cell: SaigeHeritabilityCell,
            width: saigeHeritabilityWidth,
            materialUiNoPadding: true,
            minValue: minSaigeHeritabilityValue,
            maxValue: maxSaigeHeritabilityValue,
          }
        ]
      }
      if (columnVisibilities.lambdaGc) {
        columns = [
          ...columns,
          {
            Header: "Lambda GC",
            columnGroupVisibilityAttribName: "lambdaGc",
            columns: commonPopulations.filter(pop => populationMetricsVisibilities[pop] === true).map(pop => ({
              Header: pop,
              ...getNaColumnProps(`lambda_gc_${pop}`),
              Filter: NumberRangeColumnFilter,
              filter: rangeFilterFunction,
              width: decimalNumbersColumnWidth,
            }))
          }
        ]
      }
      if (columnVisibilities.downloads) {
        columns = [
          ...columns,
          {
            Header: "Downloads",
            columnGroupVisibilityAttribName: "downloads",
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
            columnGroupVisibilityAttribName: "md5",
            columns: [
              {Header: "tsv", accessor: "md5_hex", disableFilters: true, width: columnWidth},
              {Header: "tbi", accessor: "md5_hex_tabix", disableFilters: true, width: columnWidth},
            ]
          }
        ]
      }
      return columns
    }
  , [columnVisibilities, populationMetricsVisibilities])
  const initialState = useMemo(() => {
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
    state,
    columns: outputColumns,
    setHiddenColumns,
  } = useTable<Datum>({
    columns, data: data as any,
    initialState,
    defaultColumn,
    globalFilter: fuzzyTextGlobalFilterFunction,
  },
  useFilters,
  useGlobalFilter,
  useBlockLayout)

  const headerGroupElems = headerGroups.map(headerGroup => {
    const headerElems = headerGroup.headers.map(column => {
      return (
        <TableCell align="center" {...column.getHeaderProps()}>
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

  const scrollBarSize = useMemo(() => scrollbarWidth(), [])


  const RenderRow = useCallback(({index, style}) => {
    const row = rows[index]
    prepareRow(row)
    const cellElems = row.cells.map(cell => {
      const padding = ("materialUiNoPadding" in cell.column && cell.column.materialUiNoPadding === true) ? "none" : "default"
      return (
        <TableCell {...cell.getCellProps()} padding={padding}>
          {cell.render("Cell")}
        </TableCell>
      )
    })
    return (
      <TableRow
        {...row.getRowProps({style})}
      >
        {cellElems}
      </TableRow>
    )
  }, [prepareRow, rows])



  return (
    <Layout title={`${siteConfig.title}`} description="Phenotypes">
      <header>
        <div className={`container ${phenotypesStyles.titleContainer}`}>
          <h1 className="page-title">Phenotypes</h1>
        </div>
      </header>
      <main>
        <div className={`container ${phenotypesStyles.container}`}>
          <div>
            <PhenotypeFilters
              columnVisibilities={columnVisibilities}
              setColumnVisibilities={setColumnVisibilities}
              columns={outputColumns}
              preGlobalFilteredRows={preGlobalFilteredRows}
              setGlobalFilter={setGlobalFilter}
              globalFilter={state.globalFilter}
              populationMetricsVisibilities={populationMetricsVisibilities}
              setPopulationMetricsVisibilities={setPopulationMetricsVisibilities}
            />
          </div>
          <div style={{overflowX: "auto"}}>
            <Table size="small" {...getTableProps()}>
              <TableHead>
                {headerGroupElems}
              </TableHead>
              <div {...getTableBodyProps()}>
                <FixedSizeList
                  height={700}
                  itemCount={rows.length}
                  itemSize={60}
                  width={totalColumnsWidth + scrollBarSize}
                  innerElementType={TableBody}
                >
                  {RenderRow}
                </FixedSizeList>
              </div>
            </Table>
          </div>
        </div>
      </main>
    </Layout>
  )
}

export default Phenotypes
