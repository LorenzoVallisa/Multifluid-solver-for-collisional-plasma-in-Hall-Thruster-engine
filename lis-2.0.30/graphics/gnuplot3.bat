:start

@rem Sample batch script to draw distribution of real eigenvalues using gnuplot
@rem Make executable files and run the following command:
@rem
@rem > gnuplot3.bat
@rem


@rem Set path to the utility. 

@set path="c:\Program Files (x86)\gnuplot\bin";%PATH%


@rem Run test programs.

@set bindir=..\win
@set srcdir=..\test
@%bindir%\etest5.exe %srcdir%\testmat.mtx evalues.mtx nul nul nul -e si -ie ii -ss 100


@rem Draw distribution of real eigenvalues.

@set filename=evalues.mtx
@for /f "tokens=1-5" %%a in ('findstr /v %% "%filename%"') do @(
    set size=%%a
    goto :break
    )
:break

@gnuplot.exe -e "filename='%filename%'; size=%size%" gnuplot3.plt

@del evalues.mtx

