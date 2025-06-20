function winsorize!(data, p_low=0.01, p_high=0.99)
    low = quantile(data, p_low)
    high = quantile(data, p_high)
    data[data.<low] .= low
    data[data.>high] .= high
    return data
end