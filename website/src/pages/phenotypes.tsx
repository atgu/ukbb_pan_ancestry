import React, { useCallback, useMemo, useState } from 'react'
import Layout from '@theme/Layout'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import {  useTable, useFilters, useBlockLayout  } from "react-table";
import {  SelectColumnFilter } from "../components/SelectColumnFilter";
import {  FixedSizeList } from "react-window";
import { Datum } from '../components/types';
import { scrollbarWidth } from '../components/scrollbarWidth';
import { FormGroup, FormControlLabel, Checkbox, Table, TableHead, TableRow, TableCell, TableBody, } from '@material-ui/core';
import {  DownloadLinkCell, downloadLinkCellMaxWidth } from "../components/DownloadLinkCell";
import {  CopyLinkCell, copyLinkCellMaxWidth } from "../components/CopyLinkCell";
import {  customIncludesFilterFn } from "../components/customIncludesFilterFn";
import { fuzzyTextFilterFunction } from '../components/fuzzyTextFilterFunction';
import { TextColumnFilter } from '../components/TextColumnFilter';

import data from "../data.json"
import { SliderColumnFilter } from '../components/SliderColumnFilter';
import { rangeFilterFunction } from '../components/rangeFilterFunction';
import phenotypesStyles from "./phenotypes.module.css"
import { TruncatedTextCell } from '../components/TruncatedTextCell';
import { PopulationCell } from '../components/PopulationCell';
import {  commonPopulations } from "../components/populations";



// Handle columns that can have "NA" values by replacing "NA" with `null`.
// Note: if `accessor` is a function, the `id` field must be present.
const getNaColumnProps = (fieldName) => ({
  id: fieldName,
  accessor: (row) => {
    const value = row[fieldName]
    return (value === "NA") ? null : value
  }
})

const Phenotypes = () => {
  const context = useDocusaurusContext()
  const { siteConfig = {} } = context

  const [columnVisibilities, setColumnVisibilities] = useState({
    downloads: false,
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
          disableFilters: false,
          width: 400,
          Cell: TruncatedTextCell
        },
        {
          Header: "Analysis",
          columns: [
            {
              Header: "Coding", accessor: "coding",
            },
            {
              Header: "Trait type", accessor: "trait_type",
              Filter: SelectColumnFilter,
              filter: customIncludesFilterFn,
              disableFilters: false,
            },
            {
              Header: "Phenocode", accessor: "phenocode",
              filter: fuzzyTextFilterFunction,
              Filter: TextColumnFilter,
              disableFilters: false,
            },
            {
              Header: "Sex", accessor: "pheno_sex",
              Filter: SelectColumnFilter,
              filter: customIncludesFilterFn,
              disableFilters: false,
            },
            {
              Header: "Modifier", accessor: "modifier",
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
                disableFilters: false,
                width: 400,
                Cell: TruncatedTextCell
              },
              {
                Header: "N Pops",
                accessor: "num_pops",
              },
              {
                Header: "Populations",
                accessor: "pops",
                filter: fuzzyTextFilterFunction,
                Filter: TextColumnFilter,
                disableFilters: false,
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
                Filter: SliderColumnFilter,
                filter: rangeFilterFunction,
                disableFilters: false,
              },
              {
                ...getNaColumnProps("n_cases_full_cohort_females"),
                Filter: SliderColumnFilter,
                filter: rangeFilterFunction,
                disableFilters: false,
              },
              {
                Header: "Males",
                ...getNaColumnProps("n_cases_full_cohort_males"),
                Filter: SliderColumnFilter,
                filter: rangeFilterFunction,
                disableFilters: false,
              },
              ...commonPopulations.map(pop => ({
                Header: pop,
                ...getNaColumnProps(`n_cases_${pop}`),
                Filter: SliderColumnFilter,
                filter: rangeFilterFunction,
                disableFilters: false,
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
              Filter: SliderColumnFilter,
              filter: rangeFilterFunction,
              disableFilters: false,
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
              Filter: SliderColumnFilter,
              filter: rangeFilterFunction,
              disableFilters: false,
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
              Filter: SliderColumnFilter,
              filter: rangeFilterFunction,
              disableFilters: false,
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
              },
              {
                Header: "tbi",
                accessor: "aws_link_tabix",
                Cell: DownloadLinkCell,
                maxWidth: downloadLinkCellMaxWidth,
              },
              {
                Header: "wget tsv",
                accessor: "wget",
                Cell: CopyLinkCell,
                maxWidth: copyLinkCellMaxWidth,
              },
              {
                Header: "wget tbi",
                accessor: "wget_tabix",
                Cell: CopyLinkCell,
                maxWidth: copyLinkCellMaxWidth,
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
              {Header: "tsv", accessor: "md5_hex"},
              {Header: "tbi", accessor: "md5_hex_tabix"},
            ]
          }
        ]
      }
      return columns
    }
  , [columnVisibilities])
  const initialState = useMemo(() => ({
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
      disableFilters: true,
      // Filter: DefaultColumnFilter,
    }),
    []
  )
  const {
    getTableProps, getTableBodyProps, headerGroups, rows, prepareRow,
    totalColumnsWidth,
  } = useTable<Datum>({
    columns, data: data as any,
    initialState,
    defaultColumn,
  }, useFilters, useBlockLayout)

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
