import React, { useCallback, useMemo, useState } from 'react'
import Layout from '@theme/Layout'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import {  useTable, useFilters, useBlockLayout  } from "react-table";
import {  SelectColumnFilter } from "../components/SelectColumnFilter";
import {  FixedSizeList } from "react-window";
import { Datum } from '../components/types';
import { scrollbarWidth } from '../components/scrollbarWidth';
import { FormGroup, FormControlLabel, Checkbox } from '@material-ui/core';
import {  DownloadLinkCell, downloadLinkCellMaxWidth } from "../components/DownloadLinkCell";
import {  CopyLinkCell, copyLinkCellMaxWidth } from "../components/CopyLinkCell";
import {  customIncludesFilterFn } from "../components/customIncludesFilterFn";
import { fuzzyTextFilterFunction } from '../components/fuzzyTextFilterFunction';
import { DefaultColumnFilter } from '../components/DefaultColumnFilter';

import data from "../data.json"
import { SliderColumnFilter } from '../components/SliderColumnFilter';
import { rangeFilterFunction } from '../components/rangeFilterFunction';

const pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']

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
    description: false,
    nCases: true,
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
          Header: "Analysis",
          columns: [
            {
              Header: "Coding", accessor: "coding",
            },
            {
              Header: "Phenocode", accessor: "phenocode",
              filter: fuzzyTextFilterFunction,
              Filter: DefaultColumnFilter,
              disableFilters: false,
            },
            {
              Header: "Trait type", accessor: "trait_type",
              Filter: SelectColumnFilter,
              filter: customIncludesFilterFn,
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
                Header: "Description",
                accessor: "description",
                filter: fuzzyTextFilterFunction,
                Filter: DefaultColumnFilter,
                disableFilters: false,
              },
              {
                Header: "Category",
                accessor: "category",
                filter: fuzzyTextFilterFunction,
                Filter: DefaultColumnFilter,
                disableFilters: false,
              },
              {
                Header: "N Pops",
                accessor: "num_pops",
              },
              {
                Header: "Populations",
                accessor: "pops",
                filter: fuzzyTextFilterFunction,
                Filter: DefaultColumnFilter,
                disableFilters: false,
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
              ...pops.map(pop => ({
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
            columns: pops.map(pop => ({
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
            columns: pops.map(pop => ({
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
            columns: pops.map(pop => ({
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
    // filters: [
    //   {id: "trait_type", value: undefined},
    //   {id: "n_cases_full_cohort_both_sexes", value: undefined},
    //   {id: "n_cases_full_cohort_females", value: undefined},
    //   {id: "n_cases_full_cohort_males", value: undefined},
    // ]
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
        <th {...column.getHeaderProps()}>
          {column.render("Header")}
          {filterElem}
        </th>
      )
    })
    return (
      <tr {...headerGroup.getHeaderGroupProps()}>
        {headerElems}
      </tr>
    )
  })

  const scrollBarSize = useMemo(() => scrollbarWidth(), [])


  const RenderRow = useCallback(({index, style}) => {
    const row = rows[index]
    prepareRow(row)
    const cellElems = row.cells.map(cell => {
      return (
        <div {...cell.getCellProps()}>
          {cell.render("Cell")}
        </div>
      )
    })
    return (
      <div
        {...row.getRowProps({style})}
      >
        {cellElems}
      </div>
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
        <div className="container">
          <h1 className="page-title">Results</h1>
        </div>
      </header>
      <main>
        <div className="container">
          <div className="row">
            <div className="col col--3">
              {columnFilters}
            </div>
            <div className="col col--9" style={{overflowX: "auto"}}>
              <table {...getTableProps()}>
                <thead>
                  {headerGroupElems}
                </thead>
                <div {...getTableBodyProps()}>
                  <FixedSizeList
                    height={700}
                    itemCount={rows.length}
                    itemSize={100}
                    width={totalColumnsWidth + scrollBarSize}
                  >
                    {RenderRow}
                  </FixedSizeList>
                </div>
              </table>
            </div>
          </div>
        </div>
      </main>
    </Layout>
  )
}

export default Phenotypes
