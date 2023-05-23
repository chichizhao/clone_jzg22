#!/usr/bin/env Rscript
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# x0: heterozygous SNP, x1: homozygous SNP, x2: heterozygous Indel, x3: homozygous Indel, x4: variant in all

x0 <- c(851,832,894,776,801,981,888,791,977,771,746,813,755,606,556,729,870,874,920,741,667,664,850,894,732,783,805,1046,786,904,804,717,1109,708,816,734,675,906,792,727,862,760,753,751,780,803,754,785,730,776,787,773,747,809,824,885,849,731,704,901,851,735,742,757,749,730,791,815,931,777,846,815)
x1 <- c(461,473,505,506,445,554,520,466,463,502,561,462,547,518,513,513,559,596,604,494,488,515,579,522,491,548,544,609,516,563,495,455,551,462,452,482,447,546,562,551,524,579,593,539,462,471,541,556,501,423,543,516,562,494,505,479,516,463,406,532,499,471,449,478,476,458,436,564,510,521,524,503)
x2 <- c(618,586,606,513,507,776,649,531,689,566,571,556,608,294,251,419,645,795,723,520,346,377,629,611,482,529,708,680,507,741,504,486,677,490,468,522,419,727,566,521,573,620,513,564,503,596,503,605,568,547,490,493,499,524,549,664,670,504,340,708,645,544,516,544,576,509,502,656,731,489,547,675)
x3 <- c(61,59,72,63,60,83,82,66,56,61,61,58,68,60,61,63,71,72,78,66,60,61,74,76,60,54,70,88,75,78,73,68,82,60,59,60,52,84,72,73,68,74,77,66,64,66,70,73,58,61,70,61,73,61,63,65,74,52,57,74,72,72,66,72,51,66,58,69,74,64,71,58)
x4 <- c(2016,1979,2094,1896,1861,2404,2162,1908,2225,1938,1970,1929,2030,1519,1433,1772,2167,2395,2328,1870,1621,1659,2179,2135,1805,1938,2176,2434,1929,2338,1933,1775,2419,1793,1819,1850,1682,2308,2023,1940,2071,2103,2008,1975,1838,2036,1927,2065,1886,1865,1943,1877,1925,1927,1997,2112,2160,1787,1551,2250,2102,1876,1797,1908,1909,1802,1829,2139,2264,1902,2059,2099)
list1 <- list(x0,x1,x2,x3,x4)

# Load ggplot2 library
library(ggplot2)

# Create a data frame from the list of vectors
df <- data.frame(
  Group = rep(c("Heterozygous\tSNP", "Homozygous\tSNP", "Heterozygous\tIndel", "Homozygous\tIndel", "Variant"), each = length(x0)),
  Mutations = c(x0, x1, x2, x3, x4)
)

# Create the boxplot
p <- ggplot(df, aes(x = Group, y = Mutations, fill = Group)) +
    # set the outliner color as the fill
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("#ff7f0e", "#1f77b4", "#2ca02c", "#9467bd", "#d62728")) +
  labs(
    title = NULL,
    y = "Number of mutations",
    x = NULL
  ) +
  theme(
    axis.title.x = element_text(size = rel(1.5), face = "bold"),
    axis.title.y = element_text(size = rel(4), face = "bold",vjust = 1.3),
    axis.text.x = element_text(size = rel(4), color = "black",offset(2),angle = 60, hjust = 1),
    axis.text.y = element_text(size = rel(5), color = "black"),
    axis.line = element_line(color = "black",offset(1),linewidth = 2),
    panel.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks.x = element_line(color = "black",offset(1),linewidth = 2),
    axis.ticks.y = element_line(color = "black",offset(1),linewidth = 2)
  )+   guides(fill = FALSE)

# Save the plot
path <- "/home/chichi/data/jzg_22/"
ggsave(paste0(path, "boxplot_mut.svg"), p, width = 30, height = 40, units = "cm")
# x5: heterozygous SNP ratio, x6: heterozygous Indel ratio, x7: heterozygous variant ratio
x5 <- x0/(x0+x1)
x6 <- x2/(x2+x3)
x7 <- (x0+x2)/(x0+x1+x2+x3)

list2 <- list(x5,x6,x7)

df <- data.frame(
  Group = rep(c("SNP", "Indel", "Variant"), each = length(x0)),
  Mutations = c(x5, x6, x7)
)

p <- ggplot(df, aes(x = Group, y = Mutations, fill = Group)) +
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("#ff7f0e", "#1f77b4", "#2ca02c")) +
  labs(
    title = NULL,
    y = "heterozygous rate",
    x = NULL
  ) +
  theme(
    axis.title.x = element_text(size = rel(2.5), face = "bold"),
    axis.title.y = element_text(size = rel(5), face = "bold",vjust = 1.3),
    axis.text.x = element_text(size = rel(5), color = "black",offset(2)),
    axis.text.y = element_text(size = rel(6), color = "black"),
    axis.line = element_line(color = "black",offset(1),linewidth = 2.5),
    panel.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks.x = element_line(color = "black",offset(1),linewidth = 2.5),
    axis.ticks.y = element_line(color = "black",offset(1),linewidth = 2.5)
  )+   guides(fill = FALSE)
p
ggsave(paste0(path, "boxplot_mut2.svg"), p, width = 30, height = 40, units = "cm")

# x8: heterozygous SNP, x9: homozygous SNP, x10: heterozygous Indel, x11: homozygous Indel, x12: variant in all
x8 <- c(2332,2289,2398,2177,2379,2432,2319,2265,2470,2136,2263,2240,2206,2054,1909,2141,2319,2410,2384,2313,2164,2103,2286,2374,2085,2296,2304,2651,2245,2449,2332,2228,2547,2327,2201,2283,2188,2346,2379,2397,2259,2333,2361,2172,2289,2341,2176,2320,2160,2410,2382,2280,2322,2184,2493,2427,2309,2091,2115,2442,2224,2242,2184,2246,2237,2140,2160,2279,2592,2385,2475,2421)
x9 <- c(2327,2371,2417,2296,2291,2548,2459,2270,2359,2437,2392,2357,2378,2211,2148,2313,2423,2491,2584,2348,2317,2213,2521,2441,2265,2326,2426,2512,2314,2505,2290,2199,2533,2289,2340,2301,2279,2510,2407,2398,2403,2508,2484,2320,2362,2374,2306,2411,2330,2334,2349,2302,2428,2313,2424,2431,2414,2242,2172,2486,2460,2312,2337,2336,2304,2284,2297,2481,2533,2342,2320,2475)
x10 <- c(9382,9302,9684,9152,8869,10624,9560,8624,9587,9718,9276,8954,9520,7009,7166,7853,9996,10768,10347,9253,7998,8282,9619,9437,9003,9123,9967,9175,8879,9853,8708,8886,9602,8620,9280,9239,8146,10278,8709,8497,9228,9330,9431,9174,9366,8866,8588,9060,9017,8947,8762,8991,8619,9356,8821,10298,9564,9224,7968,10111,9655,9103,9044,9335,8857,7944,8871,9910,10521,8465,9455,9848)
x11 <- c(2805,2867,2959,2864,2868,2963,2897,2734,2876,2890,2856,2859,2839,2676,2585,2834,2876,2891,2959,2807,2661,2659,2923,2868,2772,2838,2802,2953,2772,2976,2811,2756,2934,2739,2846,2894,2734,2847,2850,2863,2874,2944,2859,2828,2869,2931,2859,2940,2833,2856,2840,2813,2861,2834,2913,2849,2852,2799,2661,2944,3037,2829,2854,2868,2849,2838,2881,2975,2916,2836,2817,3013)
x12 <- x8+x9+x10+x11
list3 <- list(x8,x9,x10,x11,x12)

df <- data.frame(
  Group = rep(c("Heterozygous\tSNP", "Homozygous\tSNP", "Heterozygous\tIndel", "Homozygous\tIndel", "Variant"), each = length(x0)),
    Mutations = c(x8, x9, x10, x11, x12)
)

p <- ggplot(df, aes(x = Group, y = Mutations, fill = Group)) +
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("#ff7f0e", "#1f77b4", "#2ca02c", "#9467bd", "#d62728")) +
  labs(
    title = NULL,
    y = "heterozygous rate",
    x = NULL
  ) +
  theme(
    axis.title.x = element_text(size = rel(1.5), face = "bold"),
    axis.title.y = element_text(size = rel(4), face = "bold",vjust = 1.3),
    axis.text.x = element_text(size = rel(4), color = "black",offset(2),angle = 60, hjust = 1),
    axis.text.y = element_text(size = rel(5), color = "black"),
    axis.line = element_line(color = "black",offset(1),linewidth = 2),
    panel.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks.x = element_line(color = "black",offset(1),linewidth = 2),
    axis.ticks.y = element_line(color = "black",offset(1),linewidth = 2)
  )+   guides(fill = FALSE)
p
ggsave(paste0(path, "boxplot_mut3.svg"), p, width = 30, height = 40, units = "cm")
# x13: heterozygous SNP ratio, x14: heterozygous Indel ratio, x15: heterozygous variant ratio
x13 <- x8/(x8+x9)
x14 <- x10/(x10+x11)
x15 <- (x8+x10)/(x8+x9+x10+x11)

list4 <- list(x13,x14,x15)

df <- data.frame(
  Group = rep(c("\n SNP", "\n Indel", "\n Variant"), each = length(x0)),
    Mutations = c(x13, x14, x15)
)

p <- ggplot(df, aes(x = Group, y = Mutations, fill = Group)) +
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("#ff7f0e", "#1f77b4", "#2ca02c")) +
  labs(
    title = NULL,
    y = "heterozygous rate",
    x = NULL
  ) +
  theme(
    axis.title.x = element_text(size = rel(2.5), face = "bold"),
    axis.title.y = element_text(size = rel(5), face = "bold",vjust = 1.3),
    axis.text.x = element_text(size = rel(5), color = "black",offset(2)),
    axis.text.y = element_text(size = rel(6), color = "black"),
    axis.line = element_line(color = "black",offset(1),linewidth = 2.5),
    panel.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks.x = element_line(color = "black",offset(1),linewidth = 2.5),
    axis.ticks.y = element_line(color = "black",offset(1),linewidth = 2.5)
  )+   guides(fill = FALSE)
p
ggsave(paste0(path, "boxplot_mut4.svg"), p, width = 30, height = 40, units = "cm")