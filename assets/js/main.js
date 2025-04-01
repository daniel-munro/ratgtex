$(document).ready(function() { 
    const rn = window.genome; // Defined in the HTML pages
    const default_genome = "rn7";
    const version = "v3";

    $('#gene-search').submit(function (event) {
        var base = rn == default_genome ? "/gene#" : `/gene/${rn}#`;
        window.location.href = base + $('#search-input').val();
        return false; // prevents form submission from overriding redirect
    });

    $.ajax({
        dataType: "json",
        url: `/data/autocomplete.${version}.json`
    }).done(function (data) {
        $("#search-input").autocomplete({
            // source: data,
            source: function (request, response) {
                var results = $.ui.autocomplete.filter(data, request.term);
                response(results.slice(0, 10));
            },
            select: function(event, ui) {
                if (ui.item) {
                    $(this).val(ui.item.value);
                    $(this.form).submit();
                }
            }
        });
    });
});
