
import polars as pl

import numpy as np


def form_concentrations(experiment,printout=False):
    
    ncols = 12
    
    concentrations = [experiment['largest_conc']/(experiment['dilution']**i) for i in range(ncols-1)]
    concentrations[-1] = 0

    conc_unit = []
    units = []
    for conc in concentrations:
        cc = conc
        unit = 'ug/ml'
        if conc > 999:
            cc /= 1000
            unit = 'mg/ml'
        elif conc < 0.99:
            cc *= 1000
            unit = 'ng/ml'
        conc_unit.append(cc)
        units.append(unit)

    units[-1] = experiment['drug']

    if printout:
        print(concentrations)
        print(conc_unit)
        print(units)

    return concentrations,conc_unit,units



def read_plate(experiment,plate=0):

    df = pl.read_excel(
            experiment['file_name']+'.xlsx',
            sheet_name=f"Plate {plate+1} - Sheet1",
        )

    # Extract header row
    header = df.row(20)

    # Slice data
    df1 = df.slice(21, 151).select(df.columns[1:99])

    # Set column names
    df1.columns = header[1:99]

    # Calculate time in hours from datetime
    df1 = df1.with_columns(
        pl.col("Time").str.to_datetime().alias("Time")
    )

    t0 = df1.select(pl.col("Time").first()).item()

    df1 = df1.with_columns(
        ((pl.col("Time") - pl.lit(t0)).dt.total_seconds() / 3600).alias("t_hours")
    )

    # Calculate background
    letters = experiment['background'][0]
    numbers = experiment['background'][1]
    exclude = experiment['background'][2]

    # Build column names like A2, A3, ..., H12
    target_cols = [
        f"{letter}{number}"
        for letter in letters
        for number in numbers
        if f"{letter}{number}" not in exclude
    ]

    df2 = df1.with_columns([
        pl.concat_list(
            pl.col(target_cols).cast(pl.Float64, strict=False)
        ).list.mean().alias("backgr_mean"),

        pl.concat_list(
            pl.col(target_cols).cast(pl.Float64, strict=False)
        ).list.std().alias("backgr_std"),
    ])

    return df2


def subtract_bg_and_integrate(experiment, df, drop_first=3,time_col="t_hours"):
    
    # Form cell names
    letters = ['A','B','C','D','E','F','G','H']
    numbers = [i for i in range(1,13)]

    # Build column names like A2, A3, ..., H12
    cols = [
        f"{letter}{number}"
        for letter in letters
        for number in numbers
    ]

    # Ensure numeric
    df = df.with_columns([
        pl.col(time_col).cast(pl.Float64),
        *[pl.col(c).cast(pl.Float64, strict=False) for c in cols],
        pl.col("backgr_mean").cast(pl.Float64),
    ])

    # Background-subtracted columns
    df = df.with_columns([
        (pl.col(c) - pl.col("backgr_mean")).alias(f"{c}_bgsub")
        for c in cols
    ])

    # Cumulative trapezoidal integral over time, per column
    t = df[time_col].to_numpy()
    dt = np.diff(t)

    for c in cols:
        y = df[f"{c}_bgsub"].to_numpy()

        # Trapezoidal increments
        inc = 0.5 * (y[:-1] + y[1:]) * dt

        # Drop first noisy points
        if drop_first > 0:
            inc[:drop_first] = 0.0

        # Cumulative integral
        integ = np.concatenate([[0.0], np.cumsum(inc)])

        df = df.with_columns(
            pl.Series(f"{c}_auc", integ)
        )

    return df


def add_variant_stats(experiment, df, plate, area=False):
    """
    For each variant and each number n, compute per-row mean/std across replicate letters:
      mean( L1n, L2n, L3n ), std( ... )
    Adds columns:
      <variant>_<n>_mean, <variant>_<n>_std
    """

    numbers=range(1, 12)
    exprs = []

    for i, vname in enumerate(experiment['plates'][plate]):
        letters = experiment['variants'][i]          # FLAT list, e.g. ['F','G','H']
        letters_tag = "".join(letters)              # e.g. "FGH"

        for n in numbers:
            cols = [f"{L}{n}_auc" if area else f"{L}{n}" for L in letters]
            cols = [c for c in cols if c in df.columns]
            if not cols:
                continue

            cols_exprs = [pl.col(c).cast(pl.Float64, strict=False) for c in cols]
            k = len(cols_exprs)

            mean_expr = pl.mean_horizontal(cols_exprs)

            # sample std (ddof=1): sqrt( (sum(x^2) - sum(x)^2/k) / (k-1) )
            sum_expr = pl.sum_horizontal(cols_exprs)
            sumsq_expr = pl.sum_horizontal([e * e for e in cols_exprs])

            var_expr = pl.when(pl.lit(k) > 1).then(
                (sumsq_expr - (sum_expr * sum_expr) / k) / (k - 1)
            ).otherwise(None)

            std_expr = pl.when(var_expr.is_not_null() & (var_expr >= 0)).then(
                var_expr.sqrt()
            ).otherwise(None)

            exprs.append(mean_expr.alias(f"{vname}_{letters_tag}{n}_mean"))
            exprs.append(std_expr.alias(f"{vname}_{letters_tag}{n}_std"))


    return df.with_columns(exprs)


def get_target_times(df,targets):

    target_df = pl.DataFrame({
        "target_hour": targets
    })

    result = (
        df
        .join(target_df, how="cross")  # cross join
        .with_columns(
            (pl.col("t_hours") - pl.col("target_hour")).abs().alias("diff")
        )
        .sort("diff")
        .group_by("target_hour")
        .first()  # smallest diff per target
        .drop(["diff"])
        .sort("target_hour") 
    )

    return result

