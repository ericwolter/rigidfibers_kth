$(document).ready(function() {

    var margin = {
            top: 20,
            right: 20,
            bottom: 30,
            left: 40
        },
        width = 960 - margin.left - margin.right,
        height = 500 - margin.top - margin.bottom;

    var x = d3.scale.ordinal()
        .rangeRoundBands([0, width], .1);

    var y = d3.scale.linear()
        .rangeRound([height, 0]);

    var color = d3.scale.category10();

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select("body").append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    queue()
        .defer(d3.csv, "raw/201409081026.csv")
        .defer(d3.csv, "raw/201409081239.csv")
        .awaitAll(function(err, results) {
            console.log(results);

            color.domain(_.chain(results)
                .flatten()
                .map(function(d) {
                    return d.name
                })
                .unique()
                .value());

            var max = 0.0;
            for (var i = results.length - 1; i >= 0; i--) {
                var result = results[i];

                var total = 0.0;
                var y0 = 0.0;
                result.forEach(function(d, i) {
                    d.y0 = y0;
                    d.y1 = (y0 += +d.time);
                    total += (+d.time);
                });

                if (total > max) {
                    max = total;
                }
            };
            max *= 1.1;

            x.domain([0, results.length - 1])
            y.domain([0, max]);

            svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + height + ")")
                .call(xAxis);

            svg.append("g")
                .attr("class", "y axis")
                .call(yAxis);

            var run = svg.selectAll(".run")
                .data(results)
                .enter().append("g")
                .attr("class", "g")
                .attr("transform", function(d, i) {
                    return "translate(" + x(i) + ",0)";
                });

            run.selectAll("rect")
                .data(function(d) {
                    return d;
                })
                .enter().append("rect")
                .attr("width", x.rangeBand())
                .attr("y", function(d) {
                    return y(d.y1);
                })
                .attr("height", function(d) {
                    return y(d.y0) - y(d.y1);
                })
                .style("fill", function(d) {
                    return color(d.name);
                });

            var legend = svg.selectAll(".legend")
                .data(color.domain().slice().reverse())
                .enter().append("g")
                .attr("class", "legend")
                .attr("transform", function(d, i) {
                    return "translate(0," + i * 20 + ")";
                });

            legend.append("rect")
                .attr("x", width - 18)
                .attr("width", 18)
                .attr("height", 18)
                .style("fill", color);

            legend.append("text")
                .attr("x", width - 24)
                .attr("y", 9)
                .attr("dy", ".35em")
                .style("text-anchor", "end")
                .text(function(d) {
                    return d;
                });
        });
});