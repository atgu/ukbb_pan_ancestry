import React from "react"
import { createMuiTheme, ThemeProvider } from "@material-ui/core"
import useThemeContext from '@theme/hooks/useThemeContext';
import { PhenotypesPageContent } from "./PhenotypesPageContent";


export const PhenotypesPageContentWrapper = () => {
  const {isDarkTheme} = useThemeContext();

  const materialUiTheme = React.useMemo(
    () =>
      createMuiTheme({
        palette: {
          type: isDarkTheme ? 'dark' : 'light',
        },
        breakpoints: {
          xs: 0,
          sm: 600,
          // This is set to 966 to match the breakpoint in Infima:
          // https://github.com/facebookincubator/infima/blob/f1e03cee17347f0d54531b334fa3e4deb6a460ee/packages/core/styles/common/variables.pcss#L252
          md: 966,
          lg: 1280,
          xl: 1920,
        },
        overrides: {
          MuiTableCell: {
            root: {
              padding: "0",
            }
          }
        }
      }),
    [isDarkTheme],
  );
  return (
    <ThemeProvider theme={materialUiTheme}>
      <PhenotypesPageContent/>
    </ThemeProvider>
  )
}
