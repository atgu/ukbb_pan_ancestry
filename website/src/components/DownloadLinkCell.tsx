import React from "react"
import { Link, } from '@material-ui/core';
import {  CloudDownload as CloudDownloadIcon, } from "@material-ui/icons";

interface Props {
  value: string
}

export const DownloadLinkCell = ({value}: Props) => {
  return (
    <div style={{display: "flex", justifyContent: "center", alignItems: "center"}}>
      <Link href={value}>
        <CloudDownloadIcon/>
      </Link>
    </div>
  )
}

export const downloadLinkCellMaxWidth = 50
