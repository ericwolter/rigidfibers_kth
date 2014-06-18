$(document).ready(function() {
    vg.parse.spec('spec/ex_stacked_area.json', function(chart) { 
        chart({el:'#vis',renderer:'svg'}).update(); 
    });
});