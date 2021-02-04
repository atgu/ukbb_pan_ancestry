import React, { useCallback, useMemo, useState } from 'react'
import Layout from '@theme/Layout'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import {  useTable, useFilters, useBlockLayout, useGlobalFilter  } from "react-table";
import {  SelectColumnFilter } from "../components/SelectColumnFilter";
import {  FixedSizeList } from "react-window";
import { Datum } from '../components/types';
import { scrollbarWidth } from '../components/scrollbarWidth';
import { FormGroup, FormControlLabel, Checkbox, Table, TableHead, TableRow, TableCell, TableBody, } from '@material-ui/core';
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
import {  commonPopulations } from "../components/populations";
import { GlobalFilter } from '../components/GlobalFilter';
import { NumberRangeColumnFilter } from '../components/NumberRangeColumnFilter';


const DefaultColumnFilter = () => null

// Handle columns that can have "NA" values by replacing "NA" with `null`.
// Note: if `accessor` is a function, the `id` field must be present.
const getNaColumnProps = (fieldName) => ({
  id: fieldName,
  accessor: (row) => {
    const value = row[fieldName]
    return (value === "NA") ? null : value
  }
})
const numberRangeFilterColumnWidth = 350

const Phenotypes = () => {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context

  const [columnVisibilities, setColumnVisibilities] = useState({
    downloads: true,
    description: true,
    nCases: false,
    nControls: false,
    saigeHeritability: false,
    lambdaGc: false,
    md5: false,
  })

  const handleColumnVisibilityChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    setColumnVisibilities({
      ...columnVisibilities, [event.target.name]: event.target.checked
  }) }, [columnVisibilities])

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
          columns: [
            {
              Header: "Coding", accessor: "coding",
              width: 100,
            },
            {
              Header: "Trait type", accessor: "trait_type",
              Filter: SelectColumnFilter,
              filter: customIncludesFilterFn,
            },
            {
              Header: "Phenocode", accessor: "phenocode",
              filter: fuzzyTextFilterFunction,
              Filter: TextColumnFilter,
              Cell: TruncatedTextCell,
              width: 200,
            },
            {
              Header: "Sex", accessor: "pheno_sex",
              Filter: SelectColumnFilter,
              filter: customIncludesFilterFn,
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
            columns: [
              {
                Header: "Category",
                accessor: "category",
                filter: fuzzyTextFilterFunction,
                Filter: TextColumnFilter,
                width: 400,
                Cell: TruncatedTextCell
              },
              // {
              //   Header: "N Pops",
              //   accessor: "num_pops",
              //   disableFilters: true,
              // },
              {
                Header: "Populations",
                accessor: "pops",
                filter: fuzzyTextFilterFunction,
                Filter: TextColumnFilter,
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
            columns: [
              {
                Header: "Both sexes",
                ...getNaColumnProps("n_cases_full_cohort_both_sexes"),
                Filter: NumberRangeColumnFilter,
                filter: rangeFilterFunction,
                width: numberRangeFilterColumnWidth,
              },
              {
                Header: "Females",
                ...getNaColumnProps("n_cases_full_cohort_females"),
                Filter: NumberRangeColumnFilter,
                filter: rangeFilterFunction,
                width: numberRangeFilterColumnWidth,
              },
              {
                Header: "Males",
                ...getNaColumnProps("n_cases_full_cohort_males"),
                Filter: NumberRangeColumnFilter,
                filter: rangeFilterFunction,
                width: numberRangeFilterColumnWidth,
              },
              ...commonPopulations.map(pop => ({
                Header: pop,
                ...getNaColumnProps(`n_cases_${pop}`),
                Filter: NumberRangeColumnFilter,
                filter: rangeFilterFunction,
                width: numberRangeFilterColumnWidth,
              }))
            ]
          }
        ]
      }
      if (columnVisibilities.nControls) {
        columns = [
          ...columns,
          {
            Header: "N Controls",
            columns: commonPopulations.map(pop => ({
              Header: pop,
              ...getNaColumnProps(`n_controls_${pop}`),
              Filter: NumberRangeColumnFilter,
              filter: rangeFilterFunction,
              width: numberRangeFilterColumnWidth,
            }))
          }
        ]
      }
      if (columnVisibilities.saigeHeritability) {
        columns = [
          ...columns,
          {
            Header: "Saige heritability",
            columns: commonPopulations.map(pop => ({
              Header: pop,
              ...getNaColumnProps(`saige_heritability_${pop}`),
              Filter: NumberRangeColumnFilter,
              filter: rangeFilterFunction,
              width: numberRangeFilterColumnWidth,
            }))
          }
        ]
      }
      if (columnVisibilities.lambdaGc) {
        columns = [
          ...columns,
          {
            Header: "Lambda GC",
            columns: commonPopulations.map(pop => ({
              Header: pop,
              ...getNaColumnProps(`lambda_gc_${pop}`),
              Filter: NumberRangeColumnFilter,
              filter: rangeFilterFunction,
              width: numberRangeFilterColumnWidth,
            }))
          }
        ]
      }
      if (columnVisibilities.downloads) {
        columns = [
          ...columns,
          {
            Header: "Downloads",
            columns: [
              {
                Header: "tsv",
                accessor: "aws_link",
                Cell: DownloadLinkCell,
                maxWidth: downloadLinkCellMaxWidth,
                disableFilters: true,
              },
              {
                Header: "tbi",
                accessor: "aws_link_tabix",
                Cell: DownloadLinkCell,
                maxWidth: downloadLinkCellMaxWidth,
                disableFilters: true,
              },
              {
                Header: "wget tsv",
                accessor: "wget",
                Cell: CopyLinkCell,
                maxWidth: copyLinkCellMaxWidth,
                disableFilters: true,
              },
              {
                Header: "wget tbi",
                accessor: "wget_tabix",
                Cell: CopyLinkCell,
                maxWidth: copyLinkCellMaxWidth,
                disableFilters: true,
              },
            ]
          }
        ]
      }
      if (columnVisibilities.md5) {
        columns = [
          ...columns,
          {
            Header: "MD5",
            columns: [
              {Header: "tsv", accessor: "md5_hex", disableFilters: true},
              {Header: "tbi", accessor: "md5_hex_tabix", disableFilters: true},
            ]
          }
        ]
      }
      return columns
    }
  , [columnVisibilities])
  const initialState = useMemo(() => ({
    globalFilter: "",
    filters: [
      {id: "trait_type", value: ""},
      {id: "pheno_sex", value: ""},
      // {id: "n_cases_full_cohort_both_sexes", value: undefined},
      // {id: "n_cases_full_cohort_females", value: undefined},
      // {id: "n_cases_full_cohort_males", value: undefined},
    ]
  }), [])

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
      const filterElem = column.canFilter ? (
        <div>{column.render("Filter")}</div>
       ) : null
      return (
        <TableCell align="center" {...column.getHeaderProps()}>
          {column.render("Header")}
          {filterElem}
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
      return (
        <TableCell {...cell.getCellProps()}>
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

  const columnFilters = (
    <FormGroup row={false}>
      <FormControlLabel
        control={
          <Checkbox
            checked={columnVisibilities.downloads} name="downloads"
            onChange={handleColumnVisibilityChange}
          /> }
        label="Downloads"
      />
      <FormControlLabel
        control={
          <Checkbox
            checked={columnVisibilities.description} name="description"
            onChange={handleColumnVisibilityChange}
          /> }
        label="Description"
      />
      <FormControlLabel
        control={
          <Checkbox
            checked={columnVisibilities.nCases} name="nCases"
            onChange={handleColumnVisibilityChange}
          /> }
        label="N Cases"
      />
      <FormControlLabel
        control={
          <Checkbox
            checked={columnVisibilities.nControls} name="nControls"
            onChange={handleColumnVisibilityChange}
          /> }
        label="N Controls"
      />
      <FormControlLabel
        control={
          <Checkbox
            checked={columnVisibilities.saigeHeritability} name="saigeHeritability"
            onChange={handleColumnVisibilityChange}
          /> }
        label="Saige Heritability"
      />
      <FormControlLabel
        control={
          <Checkbox
            checked={columnVisibilities.lambdaGc} name="lambdaGc"
            onChange={handleColumnVisibilityChange}
          /> }
        label="Lambda GC"
      />
      <FormControlLabel
        control={
          <Checkbox
            checked={columnVisibilities.md5} name="md5"
            onChange={handleColumnVisibilityChange}
          /> }
        label="MD5"
      />
    </FormGroup>

  )

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
            {columnFilters}
          </div>
          <div style={{overflowX: "auto"}}>
            <Table size="small" {...getTableProps()}>
              <TableHead>
                {headerGroupElems}
                <TableRow>
                  <TableCell colSpan={visibleColumns.length}>
                    <GlobalFilter
                      preGlobalFilteredRows={preGlobalFilteredRows}
                      globalFilter={state.globalFilter}
                      setGlobalFillter={setGlobalFilter}
                    />
                  </TableCell>
                </TableRow>
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
