import { Box, Checkbox, FormControl, FormControlLabel, FormGroup, FormLabel, makeStyles, Radio, RadioGroup, Switch } from "@material-ui/core";
import React from "react"
import { commonPopulations, populationColorMapping } from "./populations";
import { FilterValue, Operator } from "./populationsFilterFunction";

const useStyles = makeStyles(() => ({
  populationCheckboxGroup: {
    display: "grid",
    gridTemplateColumns: "repeat(3, 1fr)"
  }
}))

export const PopulationsFilter = ({column}) => {
  const classes = useStyles()
  const filterValue: FilterValue | undefined = column.filterValue
  const setFilter: (value: FilterValue | undefined) => void = column.setFilter
  let content: React.ReactNode
  if (filterValue === undefined) {
    content = null
  } else {
    const handleOperatorChange = (event: React.ChangeEvent<HTMLInputElement>) => setFilter({
      ...filterValue,
      operator: (event.target as HTMLInputElement).value as Operator,
    })
    const handleSelectionChange = (event: React.ChangeEvent<HTMLInputElement>) => {
      const isChecked = event.target.checked
      const population = event.target.name
      let newPopulations: string[]
      if (isChecked === true && filterValue.populations.includes(population) === false) {
        newPopulations = [...filterValue.populations, population]
      } else if (isChecked === false && filterValue.populations.includes(population) === true) {
        newPopulations = filterValue.populations.filter(pop => pop !== population)
      } else {
        newPopulations = filterValue.populations
      }
      setFilter({
        ...filterValue,
        populations: newPopulations,
      })
    }

    const checkboxes = commonPopulations.map(population => {
      const color = populationColorMapping.get(population)!
      return (
        <FormControlLabel key={population}
          label={population}
          style={{color}}
          control={
            <Checkbox
              style={{color}}
              checked={filterValue.populations.includes(population)}
              onChange={handleSelectionChange}
              name={population}
            />
          }
        />
      )
    })
    content = (
      <>
        <FormControl component="fieldset">
          <FormLabel component="legend">Filter Operator</FormLabel>
          <RadioGroup value={filterValue.operator} onChange={handleOperatorChange} row={true}>
            <FormControlLabel value={Operator.AND} control={<Radio/>} label="AND"/>
            <FormControlLabel value={Operator.OR} control={<Radio/>} label="OR"/>
          </RadioGroup>
        </FormControl>
        <FormControl component="fieldset">
          <FormLabel component="legend">Selection</FormLabel>
          <FormGroup className={classes.populationCheckboxGroup}>
            {checkboxes}
          </FormGroup>
        </FormControl>
      </>
    )
  }
  const handleSwitchChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const isChecked = event.target.checked
    if (isChecked) {
      const populationFilterValue: FilterValue = {
        operator: Operator.AND,
        populations: commonPopulations,
      }
      setFilter(populationFilterValue)
    } else {
      setFilter(undefined)
    }
  }
  return (
    <Box display="flex" m={1}>
      <FormControl component="fieldset">
        <FormControlLabel
          control={
            <Switch
              checked={filterValue !== undefined}
              onChange={handleSwitchChange}
            />
          }
          label="Combination Selection"
        />
        {content}
      </FormControl>
    </Box>
  )
}
