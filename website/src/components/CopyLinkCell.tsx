import React from "react"
import { IconButton, makeStyles, } from '@material-ui/core';
import {  FileCopy as FileCopyIcon } from "@material-ui/icons";

const useStyles = makeStyles({
  root: {
    display: "flex",
    justifyContent: "center",
    alignItems: "center",
    height: "100%",
  }
})

interface Props {
  value: string
}

export const CopyLinkCell = ({value}: Props) => {
  const classes = useStyles()
  const onClick = () => {
    // TODO: actually copy to clipboard here
  }
  return (
    <div className={classes.root}>
      <IconButton color="primary" onClick={onClick}>
        <FileCopyIcon/>
      </IconButton>
    </div>
  )
}

export const copyLinkCellMaxWidth = 80
