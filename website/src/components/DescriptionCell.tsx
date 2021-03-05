import { Link, makeStyles, Theme, Tooltip } from "@material-ui/core";
import React from "react"
import clsx from "clsx"
import {  height } from "./BarChartCell";

export const width = 300;

const useStyles = makeStyles((theme: Theme) => ({
  root: {
    display: "grid",
    gridTemplateRows: "repeat(3, 1fr)",
    gridTemplateColumns: "1fr 1fr 2fr",
    alignItems: "center",
    alignContent: "center",
    height: `${height}px`,
    padding: theme.spacing(0.5),
  },
  // https://css-tricks.com/snippets/css/truncate-string-with-ellipsis/
  truncatedTextRight: {
    width: `${width - 2 * theme.spacing(0.5)}px`,
    whiteSpace: "nowrap",
    overflow: "hidden",
    textOverflow: "ellipsis",
  },
  // https://stackoverflow.com/a/9793669/7075699
  truncatedTextLeft: {
    whiteSpace: "nowrap",
    overflow: "hidden",
    textOverflow: "ellipsis",
    width: `${width - 2 * theme.spacing(0.5)}px`,
    direction: "rtl",
    textAlign: "left",
  },
  name: {
    gridColumn: "1 / 4",
    gridRow: "1",
  },
  category: {
    gridColumn: "1 / 4",
    gridRow: "2",
    fontSize: "0.75rem",
  },
  small: {
    fontSize: "0.7rem",
    textAlign: "center"
  },
  sex: {
    gridColumn: "1",
    gridRow: "3",
  },
  code: {
    gridColumn: "2",
    gridRow: "3",
    display: "block",
  },
  traitType: {
    gridColumn: "3",
    gridRow: "3",
  }
}))

export interface Description {
  name: string
  category: string
  sex: string | undefined
  code: string | undefined
  traitType: string | undefined
}

interface Props {
  value: Description
}

export const DescriptionCell = ({value}: Props) => {
  const classes = useStyles()
  const {name, category, sex, code, traitType} = value;
  const sexElem = sex ? (
    <div className={clsx(classes.sex, classes.small)}>
      Sex: {sex === "both_sexes" ? "both" : sex}
    </div>
  ) : null
  const codeElem = code ? (
    <div className={clsx(classes.code, classes.small)} >
      Code:{" "}
      <Link
        href={`https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=${code}`}
        target="_blank" rel="noopener noreferrer"
      >
        {code}
      </Link>
    </div>
  ) : null
  const traitTypeElem = traitType ? (
    <div className={clsx(classes.traitType, classes.small)}>
      Trait type: {traitType}
    </div>
  ) : null
  return (
    <div className={classes.root}>
      <Tooltip title={name} placement="top-start">
        <div className={clsx(classes.truncatedTextRight, classes.name)}>{name}</div>
      </Tooltip>
      <Tooltip title={category} placement="top">
        <div className={clsx(classes.truncatedTextLeft, classes.category)}>{category}</div>
      </Tooltip>
      {sexElem}
      {codeElem}
      {traitTypeElem}
    </div>
  )
}
