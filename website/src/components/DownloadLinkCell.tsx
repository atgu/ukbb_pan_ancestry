import React from "react"
import { Link, makeStyles, } from '@material-ui/core';
import {  CloudDownload as CloudDownloadIcon, } from "@material-ui/icons";

const useStyles = makeStyles({
  root: {
    display: "flex",
    justifyContent: "center",
    alignItems: "center",
    height: "100%"
  }
})

interface Props {
  value: string
}

export const DownloadLinkCell = ({value}: Props) => {
  const classes = useStyles()
  return (
    <div className={classes.root}>
      <Link href={value}>
        <CloudDownloadIcon/>
      </Link>
    </div>
  )
}

export const downloadLinkCellMaxWidth = 70
