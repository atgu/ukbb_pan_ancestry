import React from "react"
import {  PopulationColorLegend } from "./PopulationColorLegend";

export const PopuplationHeader = ({column}) => {
  return (
    <div {...column.getHeaderProps()}>
      Populations
      <PopulationColorLegend iconFontSize="1rem" labelFontSize="0.75rem"/>
    </div>
  )
}
